import numpy as np
from pymatgen.transformations.advanced_transformations  import CubicSupercellTransformation
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator


"""
 DEPENDENCIES:
 impurity generator now in a  pymatgen extension
 (i) pymatgen-analysis-defects (pymatgen>=2022.10.22)
 (ii) numpy
"""

class SCgenerators:
    """
    Generates a nearly cubic supercell (SC) for convergence checks.
    Inserts a muon in the supercell at a Voronoi interstitial site.
    One of it methods initializes the supercell generation and the other 
    re-initializes generation of a larger supercell-size than the former.
    
    Param: 
        py_struc: A pymatgen "unitcell" structure instance 
                  This is used to create the supercell.  
    """
    
    @staticmethod
    def gen_nearcubic_SC(
        py_struc,
        min_atoms,
        max_atoms,
        min_length
    ):
        """
        Function that generates the nearly cubic supercell (SC).
        
        Params:
            py_struc         : The pymatgen structure instance
            min_atoms        : Integer-->Min number of atoms allowed in SC
            max_atoms        : Integer-->Max number of atoms allowed in SC 
            min_length       : Integer-->Min length of the smallest SC lattice vector 
         
        Returns: 
            A nearly cubic SC structure and an array of the SC grid size  
        """
        
        CST = CubicSupercellTransformation(
            min_atoms  = min_atoms,
            max_atoms  = max_atoms,
            min_length = min_length,
            force_diagonal = True
        ) 
        
        py_SCstruc = CST.apply_transformation(py_struc)
        SC_matrix    = np.dot(py_SCstruc.lattice.matrix,py_struc.lattice.inv_matrix)
        #SC_matrix   = np.array(SC_matrix).astype(int)
        SC_matrix    = np.asarray(np.round_(SC_matrix), dtype = 'int')
        
        return py_SCstruc, SC_matrix
    
    
    """PS: func below not used within the class and can be removed if not necessary in future!"""
    @staticmethod
    def gen_SC_from_grid(py_struc, SC_matrix):
        """
        Function that generates supercell structure for a given grid size.
        """
        py_SCstruc = py_struc.copy()
        py_SCstruc.make_supercell(SC_matrix)
        return py_SCstruc

    
    
    @staticmethod
    def append_muon_to_SC(py_SCstruc,SC_matrix,mu_frac_coord):
        """
        Add the muon as a hydrogen atom to the supercell (SC).
        
        Params:
            py_SCstruc    : The pymatgen supercell structure
            SC_matrix           : array-->the SC grid size
            mu_frac_coord     : array-->Interstitial site scaled in units 
        
        Returns: 
            A Pymatgen supercell structure that has the muon as a H atom at a Voronoi site
            
        
        """
        
        mu_frac_coord_SC  = np.dot(mu_frac_coord,np.linalg.inv(SC_matrix))
        py_SCstruc_withmu = py_SCstruc.copy()
        
        """ what if a H specie is in the structure object? """
        try:
            py_SCstruc_withmu.append(
                species = "H",
                coords  = mu_frac_coord_SC, 
                coords_are_cartesian = False, 
                validate_proximity   = True
            )
        except ValueError:
            raise SystemExit(
                'ValueError:The muon is too close to an existing site!, change muon site. Exiting....'
            ) from None
            
        return py_SCstruc_withmu


    def __init__(self,py_struc):
        
        self.py_struc       = py_struc
        #self.py_SCstruc     = None
        #self.mu_frac_coord  = None
    
    
    def initialize(self):
        """
        This func initializes the first supercell (SC) generation
        with the following conditions;
        
        min_atoms  : number of atoms in the given struc + 1 
        max_atoms  : number of atoms in the given struc*(2**3)
                    This limits the SC generation to 8 times of the given cell.            
        min_length : Min length of the smallest SC lattice vector + 1
        
        Returns: 
            A Pymatgen supercell structure that has the muon as a H atom at a Voronoi site
        """
        
        min_atoms  = self.py_struc.num_sites+1   
        max_atoms  = self.py_struc.num_sites*(2**3)
        min_length = np.min(self.py_struc.lattice.abc)+1

        py_SCstruc,SC_matrix = self.gen_nearcubic_SC(
            self.py_struc,
            min_atoms,
            max_atoms,
            min_length
        )
        
        
        """ This check will be out after testing as max_atom_num is not user defined"""
        if min_atoms-1 >= py_SCstruc.num_sites:
            raise Exception(
                'Supercell not created: Revisit the SC max_atom_num and min_length conditions.'
            ) 
        
        """ get a Voronoi interstitial site for the muon impurity, CALL NICHE?  """
        vig = VoronoiInterstitialGenerator()
        mu_frac_coord = list(vig._get_candidate_sites(self.py_struc))[0][0]    
        
        py_SCstruc_with_mu = self.append_muon_to_SC(
            py_SCstruc,
            SC_matrix,
            mu_frac_coord)
        
        return py_SCstruc_with_mu,SC_matrix,mu_frac_coord
    
    def re_initialize(self,py_SCstruc_with_mu, mu_frac_coord, iter_num): 
        """
        This function re-initializes the generation of a larger supercell-size in a loop 
        when a condition is not met after the first initialization above.
        
        Param:
            iter_num : Integer--> iteration number in the loop
        
        Returns: 
            A Pymatgen supercell structure that has the muon as a H atom at a Voronoi site
        """
        
        min_atoms  = py_SCstruc_with_mu.num_sites
        max_atoms  = self.py_struc.num_sites*((2+iter_num)**3)
        min_length = np.min(py_SCstruc_with_mu.lattice.abc)+1
        
        
        """ This check will be out after testing"""
        if min_atoms > max_atoms:
            raise ValueError('min_atoms > max_atom while re-generating supercell, check (restart_num).')
        
        py_SCstruc,SC_matrix = self.gen_nearcubic_SC(
            self.py_struc,
            min_atoms,
            max_atoms,
            min_length
        )

        """ This check will be out after testing"""
        if min_atoms-1 >= py_SCstruc.num_sites:
            raise Exception('Supercell not created: Revisit the SC max_atom_num and min_length conditions.') 
        

        py_SCstruc_with_mu = self.append_muon_to_SC(
            py_SCstruc,
            SC_matrix,
            mu_frac_coord)
        
        return py_SCstruc_with_mu,SC_matrix




import argparse
from pymatgen.core import Structure



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate nearly cubic supercell")
    parser.add_argument(
        "--iter_num",
        metavar="N",
        type=int,
        default=1,
        help="iteration number in the supercell convergence loop",
    )
    parser.add_argument("input_structure")
    
    args = parser.parse_args()
    
    iter_num = args.iter_num
    
    # load structure with pymatgen
    py_struc = Structure.from_file(args.input_structure)
    
    sg    = SCgenerators(py_struc)
    
    #initialize the caluclations
    py_SCstruc_mu2,SC_matrix,mu_frac_coord=sg.initialize()
    py_SCstruc_mu2.to(filename="positions.cif".format())
    
    # while and if loop then depending on workchain usage
    #py_SCstruc_mu2,SC_matrix=sg.re_initialize(py_SCstruc_mu2,mu_frac_coord,iter_num)
