#ALOT OF WORK STILL HAVE TO BE DONE 
# include pseudos as input, exclude from the PWworkcahin or find a way around to inlcude muon site not present in the initial input structure
#COMMENTS
#CLEAN-UP
#RE-LOGICING DUPLICATED CODES
#MAKE PARSER, PLUGIN
#RUN MORE EXAMPLES
#RESTRUCTURE
#BETTER VARIABLE AND CLASSES NAME
#BUT THEN THE WORKCHAIN STRUCTURE/LOGIC WILL BE MUCH SIMILAR
#----------------------------
from aiida import orm
from aiida.engine import ToContext, WorkChain, calcfunction,  workfunction
from aiida.orm import AbstractCode, Int, Code, Str, Float, Bool, Group, List
from aiida.orm import Dict, KpointsData, StructureData, ArrayData,load_code, load_group,load_node
#from aiida.plugins.factories import CalculationFactory
from aiida.plugins import CalculationFactory, WorkflowFactory, DataFactory
from aiida import load_profile
from aiida.engine import run, submit, CalcJob, launch
from aiida_quantumespresso.utils.mapping import prepare_process_inputs, update_mapping
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.nodes.data.upf import get_pseudos_from_structure
from aiida_pseudo.data.pseudo.upf import UpfData
import numpy as np
from aiida.engine import if_, while_, return_
from supcgen import SCgenerators
from chkconv import check_SC_convergence
load_profile()


@calcfunction
def gensup(aiida_structure):
    py_struc=aiida_structure.get_pymatgen_structure()
    py_SCstruc = py_struc.copy()
    py_SCstruc.make_supercell([2,2,2])
    py_SCstruc_mu = py_SCstruc.copy()
    py_SCstruc_mu.append(species="H", coords=[0.125,0.125,0.125], coords_are_cartesian= False, validate_proximity=True)
    StructureData = DataFactory('structure')
    aiida_structure2 = StructureData(pymatgen = py_SCstruc_mu)
    return aiida_structure2



@calcfunction
def init_supcgen(aiida_struc):
    py_struc=aiida_struc.get_pymatgen_structure()
    #
    scg=SCgenerators(py_struc)
    py_SCstruc_mu=scg.initialize()
    #
    StructureData = DataFactory('structure')
    aiida_SCstruc = StructureData(pymatgen = py_SCstruc_mu)
    
    return aiida_SCstruc
    
    


@calcfunction
def re_init_supcgen(aiida_struc,aiida_SCstruc,iter_num):
    py_struc=aiida_struc.get_pymatgen_structure()
    py_SCstruc=aiida_SCstruc.get_pymatgen_structure()
    
    #
    scg=SCgenerators(py_struc)
    py_SCstruc_mu=scg.re_initialize(py_SCstruc,mu_frac_coord,iter_num.value)
    #
    StructureData = DataFactory('structure')
    aiida_SCstructure = StructureData(pymatgen = py_SCstruc_mu)
    return aiida_SCstructure

#@workfunction
@calcfunction
def conv_achieved(aiida_structure,traj_out):
    atm_forc=traj_out.get_array('forces')[0]
    atm_forces=np.array(atm_forc)
    ase_struc = aiida_structure.get_ase()
    #
    csc=check_SC_convergence(ase_struc,atm_forces)
    cond=csc.apply_first_crit()
    cond2=csc.apply_2nd_crit()
    #
    BoolData = DataFactory('core.bool')
    if cond==True and all(cond2):
        return BoolData(True)
    else:
        return BoolData(False)


PwCalculation = CalculationFactory('quantumespresso.pw')


class muSConvWorkChain(WorkChain):
    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.input("structure", valid_type=orm.StructureData,required=True, help='Input initial structure')
        spec.expose_inputs(PwCalculation, namespace='pwscf',exclude=('structure'))   #use the  pw calcjob
        
        spec.outline(cls.init_supcell_gen,
                     cls.run_pw_scf,
                     cls.inspect_run_get_forces,
                     while_(cls.continue_iter)(
                         cls.increment_n_by_one,
                         if_(cls.iteration_num_not_exceeded)(
                             cls.get_larger_cell,
                             cls.run_pw_scf,
                             cls.inspect_run_get_forces
                         )
                         .else_(
                             cls.exit_max_iteration_exceeded,
                         )
                     ),
                     cls.set_outputs,
                    )
        
        

        spec.output('Converged_supercell', valid_type=StructureData)
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_SCF',message='one of the PwCalculation subprocesses failed')
        spec.exit_code(702, 'ERROR_NUM_CONVERGENCE_ITER_EXCEEDED',message='Max number of supercell convergence reached ')
    
    def init_supcell_gen(self):
        self.ctx.n = 0
        self.ctx.max_it_n = 2 #decide in meeting
        #self.ctx.sup_struc_mu,self.ctx.vor_site = init_supcgen(self.inputs.structure)
        self.ctx.sup_struc_mu = init_supcgen(self.inputs.structure)
    
    def run_pw_scf(self):
        #
        inputs = AttributeDict(self.exposed_inputs(PwCalculation, namespace='pwscf'))
        inputs.structure=self.ctx.sup_struc_mu
        
        running = self.submit(PwCalculation, **inputs)
        self.report(f'running SCF calculation {running.pk}')
        #self.to_context(calculation_run=running)
        
        return ToContext(calculation_run=running)
    
    def inspect_run_get_forces(self):
        calculation = self.ctx.calculation_run
        if not calculation.is_finished_ok:
            self.report('PwCalculation<{}> failed with exit status {}' .format(calculation.pk, calculation.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF
        else:
            
            self.ctx.traj_out=calculation.outputs.output_trajectory
            
    def continue_iter(self):
        #check convergence and decide if to continue the loop
        conv_res=conv_achieved(self.ctx.sup_struc_mu,self.ctx.traj_out)
        return conv_res.value == False    
        ##This implies
        #if conv_res.value == False:
        #    return True
        #else:
        #    return False   
    
    def increment_n_by_one(self):
        self.ctx.n += 1
        
    def iteration_num_not_exceeded(self): 
        return self.ctx.n <= self.ctx.max_it_n
    
    
    def get_larger_cell(self):
        self.ctx.sup_struc_mu = re_init_supcgen(
            self.inputs.structure,
            self.ctx.sup_struc_mu,
            self.ctx.n
        )
        
    
    def exit_max_iteration_exceeded(self):
        self.report('Exiting muSConvWorkChain, Coverged supercell NOT achieved, next iter num <{}> is greater than max iteration number {}' .format(self.ctx.n, self.ctx.max_it_n))
        return self.exit_codes.ERROR_NUM_CONVERGENCE_ITER_EXCEEDED
        
    
    def set_outputs(self):
        self.report('Setting Outputs')
        self.out('Converged_supercell',self.ctx.sup_struc_mu)
# END








from pymatgen.io.cif import CifParser
if __name__ == '__main__':
    parser = CifParser("Si.cif")
    py_struc = parser.get_structures()[0]
    StructureData = DataFactory('structure')
    aiida_structure = StructureData(pymatgen = py_struc)


    builder=muSConvWorkChain.get_builder()
    structure = aiida_structure
    builder.structure = structure
    codename = 'pw7_0@localhost_serial'
    code = Code.get_from_string(codename)
    #builder = code.get_builder()
    builder.pwscf.code=code
    #
    family = load_group('SSSP/1.2/PBE/efficiency')
    pseudos = family.get_pseudos(elements=('Si', 'H'))
    #pseudos = family.get_pseudos(structure=structure)
    builder.pwscf.pseudos = pseudos

    Dict = DataFactory('dict')
    parameters = {
    'CONTROL': {
    'calculation': 'scf',
    'restart_mode': 'from_scratch',
    'tstress':True,
    'tprnfor':True,
    },
    'SYSTEM': {
    'ecutwfc': 30.,
    'ecutrho': 240.,
    'tot_charge': 1.0,
    #'nspin': 2,
    'occupations':'smearing',
    'smearing':'cold',
    'degauss':0.01,

    },
    'ELECTRONS': {
    'conv_thr': 1.e-6,
    'electron_maxstep':300,
    'mixing_beta':0.3,
    }
    }

    builder.pwscf.parameters = Dict(dict=parameters)
    #
    KpointsData = DataFactory('array.kpoints')
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([1,1,1])
    builder.pwscf.kpoints = kpoints
    settings_dict={}
    settings_dict={'gamma_only': True}
    builder.pwscf.settings = Dict(dict=settings_dict)
    builder.pwscf.metadata.description = 'a PWscf  test SCF'
    builder.pwscf.metadata.options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine' : 1}
    #
    results, node = run.get_node(builder)
    node.exit_status #to check if the calculation was successful
    print(results) # from #results, node = run.get_node(builder)
