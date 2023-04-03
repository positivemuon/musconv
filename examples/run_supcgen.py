from musConv.supcgen import SCgenerators

from pymatgen.core import Structure

if __name__ == "__main__":
	# load structure with pymatgen
    py_struc = Structure.from_file("LiF.cif")
    
    sg    = SCgenerators(py_struc)
    
    #initialize the caluclations
    #py_SCstruc_mu2,SC_matrix,mu_frac_coord=sg.initialize(min_length)
    py_SCstruc_mu2,SC_matrix,mu_frac_coord=sg.initialize()
    py_SCstruc_mu2.to(filename="positions.cif".format())
    #print(SC_matrix)
    
    # while and if loop then depending on workchain usage
    #py_SCstruc_mu2,SC_matrix=sg.re_initialize(py_SCstruc_mu2,mu_frac_coord)
