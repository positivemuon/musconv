from musConv.chkconv import chkSCconvergence
from ase.io import read
import numpy as np

if __name__ == "__main__":

	# load structure with ase and then forces
    ase_struc = read("LiF_p1.cif")
    atf = np.loadtxt("LiF_p1_forces.txt")
    
    #call the func
    csc   = chkSCconvergence(ase_struc,atf)
    cond  = csc.apply_first_crit()
    cond2 = csc.apply_2nd_crit()
    print(f"Convergence of 1st criteria is {cond}, while 2nd criteria is {cond2}")
