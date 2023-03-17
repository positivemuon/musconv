from musConv.aiida_muSConvWorkChain import muSConvWorkChain
from aiida import orm
from aiida.plugins import  DataFactory
from aiida import load_profile
from aiida.engine import run, submit, ToContext
load_profile()

from pymatgen.io.cif import CifParser


if __name__ == '__main__':
    parser = CifParser("Si.cif")
    py_struc = parser.get_structures()[0]
    aiida_structure = orm.StructureData(pymatgen = py_struc)


    builder=muSConvWorkChain.get_builder()
    structure = aiida_structure
    builder.structure = structure
    builder.kpoints_distance = orm.Float(0.401)
    builder.pseudofamily = orm.Str('SSSP/1.2/PBE/efficiency')

    """N:B the pseudos and kpoints are no longer inputs in pwworkchain, 
       already taken care of in the workchain
    """

    codename = 'pw7_0@localhost_serial' #edit pw code name
    code = orm.Code.get_from_string(codename)
    builder.pwscf.code = code

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
    builder.pwscf.metadata.description = 'a PWscf  test SCF'
    builder.pwscf.metadata.options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine' : 1}
    #
    results, node = run.get_node(builder)
    
    """
    #or
    #node = submit(builder)
    #node.exit_status #to check if the calculation was successful

    # Get  converged supercell results with run
    #print(results) # from #results, node = run.get_node(builder)
    #py_conv_struct=results['Converged_supercell'].get_pymatgen_structure()
    #py_conv_struct.to(filename="supercell_withmu.cif".format())
    #Sc_matrix=results['Converged_SCmatrix'].get_array('SC_matrix')

    #get results  with submit
    #res=orm.load_node(node.pk)
    #py_conv_struct=res.outputs['Converged_supercell'].get_pymatgen_structure()
    #py_conv_struct.to(filename="supercell_withmu.cif".format())
    #Sc_matrix=res.outputs['Converged_SCmatrix'].get_array('SC_matrix')
    """
