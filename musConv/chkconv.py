import numpy as np
from ase.units import Bohr,Rydberg
from scipy.optimize import curve_fit


class chkSCconvergence:
    """
    Checks if a supercell (SC) size is converged for muon site calculations
    using results of atomic forces from a one shot SCF calculation.
    """
    
    @staticmethod
    def exp_fnc(xdata,A,B):
        """
        An exponential decay function with; 
        """
        return A*np.exp(-B*xdata)
    
    
    @staticmethod
    def min_SCconv_dist(y_C,A,B):    
        """ 
        Inverse of the exp func
        """
        return np.log(y_C/A)/(-B)

    
    
    def __init__(
        self,
        ase_struc,
        atomic_forces,
        mu_num_spec:  int or str   = 1 or 'H',
        conv_thr:     float = 1e-03*Rydberg/Bohr, #from au to eV/A
        max_force_thr:float = 0.06*Rydberg/Bohr,   #from au to eV/A
        mnasf        : int  = 4
    ):
        
        """
        Params:
        ase_struc  : An ASE Atom structure object
        
        atomic_forces  : ndarray--> array of atomic_forces in eV/A
                         Length and order of array is same as the atomic positions                      
        mu_num_spec    : Integer or String --> atomic number or Specie symbol for the muon
                         Default: 1 or 'H'                         
        conv_thr       : Float --> Converged Force threshold. 
                         Default = 1e-03*Rydberg/Bohr #from au to eV/A                         
        max_force_thr  : Float --> Max force considered in fitting
                         Default = 0.06*Rydberg/Bohr #from au to eV/A
        mnasf          : Int   -->  Minimum number of atoms sufficient for fit               
        """
        
        self.ase_struc       = ase_struc
        self.atomic_forces   = atomic_forces 
        self.mu_num_spec     = mu_num_spec 
        self.conv_thr        = conv_thr 
        self.max_force_thr   =  max_force_thr
        self.mnasf           =  mnasf
        
        assert(len(self.ase_struc.numbers)==len(self.atomic_forces)),'No. of atoms not equal to number of forces' 
        
        
        """ Check and get muon index"""
        if isinstance(self.mu_num_spec, int) and self.mu_num_spec in self.ase_struc.numbers:
            mu_idd = [atom.index for atom in self.ase_struc if atom.number == self.mu_num_spec]
            
        elif isinstance(self.mu_num_spec, str) and  self.mu_num_spec in set(self.ase_struc.get_chemical_symbols()):
            mu_idd = [atom.index for atom in self.ase_struc if atom.symbol == self.mu_num_spec]
            
        else:
            raise ValueError(f"{mu_num_spec} is not in the specie or atomic number list")  
              
        if len(mu_idd) > 1:
            raise Exception('Provided muon specie or atomic number has more than one muon in the structure')
        self.mu_id = mu_idd[0]
        
        
        """magnitude of atomic forces"""
        self.atm_forces_mag = [np.sqrt(x.dot(x)) for x in self.atomic_forces]  
     
    
    def apply_first_crit(self): 
        """
        Implements the first convergence criteria;
        Convergence is achieved if one of the forces in the supercell 
        (SC) is less than the force convergence threshold
        """
        
        # remove forces on muon  and 0. forces at boundary due to symmetry if present
        no_mu_atm_forces_mag1 = [v for i,v in enumerate(self.atm_forces_mag) if i != self.mu_id and v != 0.] 
        
        if min(no_mu_atm_forces_mag1) < self.conv_thr: 
            print("First SC convergence Criteria achieved")
            return True
        else:         
            print("First SC convergence Criteria  NOT achieved")
            return False
        
    def apply_2nd_crit(self): 
        """
        Implements the second convergence criteria:
        Forces for each atomic specie are fitted to an exponential 
        decay of their respective atomic position distance from the muon.
        Convergence is achieved when the max relative distance is less than
        the minimum relative distance obtained fro the fit parameters.
        
        """    
        atm_indxes = [atom.index for atom in self.ase_struc] 
        
        atm_dist   = self.ase_struc.get_distances(self.mu_id, atm_indxes, mic = True, vector = False) 
        
        specie_num = len(set(self.ase_struc.numbers))  
        specie_set = set(self.ase_struc.get_chemical_symbols()) 
        mu_symb    = self.ase_struc.symbols[self.mu_id]    
        
        # Remove muon specie from specie set 
        specie_set = [i for i in specie_set if i !=mu_symb] 
        
        cond = []
        for i in range(0,specie_num-1):
            specie_index = [atom.index for atom in self.ase_struc if atom.symbol == specie_set[i]]
            
            if len(specie_index) >= self.mnasf:
                specie_dist  = [atm_dist[i] for i in specie_index if self.atm_forces_mag[i] < self.max_force_thr and self.atm_forces_mag[i] != 0.] 
                specie_force = [self.atm_forces_mag[i] for i in specie_index if self.atm_forces_mag[i] < self.max_force_thr and self.atm_forces_mag[i] != 0.]

                try:
                    par,cov  = curve_fit(self.exp_fnc,specie_dist, specie_force)
                except (ValueError, RuntimeError) as err:           
                    #raise Exception("Check force data, maybe the data does not decay exponentially")
                    print("Check force data, maybe the data does not decay exponentially with:", err)
                    cond.append(False)
                    continue

                
                """fit  and data check, better conditions?"""
                stder= np.sqrt(np.diag(cov))   
                if stder[0] > par[0] or stder[1] > par[1]:
                    print(f"Check force data and fit on specie {specie_set[i]}, maybe does not decay exponentially")
                    cond.append(False)
                    continue
                    
                
                """find min distance  required for convergence"""
                min_conv_dist = self.min_SCconv_dist(self.conv_thr,par[0],par[1]) 
                #print(f"Max mu-atom distance in the SC is {max(specie_dist)}, Min distance required for conv is {min_conv_dist}")
                
                
                """TO DO: print lines to be removed"""
                if max(specie_dist) >= min_conv_dist:
                    print('Second SC convergence Criteria achieved on specie--> {}'.format(specie_set[i]))
                    cond.append(True)
                else:
                    print(f"For specie {specie_set[i]} the 2nd SC convergence is NOT achieved, min dist required is {min_conv_dist} Ang ")
                    cond.append(False)
                    
            else:
                print('The current SC size is NOT sufficient for convergence checks on specie--> {}'.format(specie_set[i]))
        
        if not cond:
            cond.append(False)
                
        return cond



import argparse
from ase.io import read

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=" Check supercell convergence")
    
    parser.add_argument(
        "--mu_num_spec",
        metavar="A",
        type=int or str,
        default=1 or "H",
        help="Muon atomic number of specie index",
     )
    
    parser.add_argument(
        "--conv_thr",
        metavar="C",
        type=float,
        default=1e-03*Rydberg/Bohr,
        help="Converged Force threshold",
    )
    
    parser.add_argument(
        "--max_force_thr",
        metavar="F",
        type=float,
        default=0.06*Rydberg/Bohr,
        help="Max force considered in fitting",
    )
           
    parser.add_argument("input_structure")
    parser.add_argument("atomic_forces")
    
    args = parser.parse_args()
    
    mu_num_spec   = args.mu_num_spec
    conv_thr      = args.conv_thr
    max_force_thr = args.max_force_thr
    
    # load structure with ase and then forces
    ase_struc = read(args.input_structure)
    atf=np.loadtxt(args.atomic_forces)
    
    #call the func
    csc   = chkSCconvergence(ase_struc,atf)
    cond  = csc.apply_first_crit()
    cond2 = csc.apply_2nd_crit()
    print(f"Convergence of 1st criteria is {cond}, while 2nd criteria is {cond2}")
