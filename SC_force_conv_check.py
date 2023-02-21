#
import numpy as np
from ase.io import read
from ase.units import Bohr,Rydberg
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#--------------------------CLASS 3: CHECK SUPERCELL CONVERGENCE BY ANALYZING FORCES IN PWOUTPUT FILE FROM AIIDA------------

def exp_fnc(xdata,A,B):
    return A*np.exp(-B*xdata)

def min_SCconv_dist(conv_thr,A,B):    #find min xdata for conv
    return np.log(conv_thr/A)/(-B)

def fit_curve(e_fnc,xdata,ydata):
    parameters, covariance = curve_fit(e_fnc, xdata, ydata)
    return parameters, covariance

def plot_forces_vs_dist_frm_mu(xdata,ydata,fit_par,sp_labl,color):
    xdata_sorted=np.sort(xdata)
    plt.plot(xdata_sorted, exp_fnc(xdata_sorted,fit_par[0],fit_par[1]), linestyle='dashed',color=color)
    plt.legend([sp_labl],prop={'size':15})
    plt.plot(xdata, ydata,color=color,linestyle = 'None',marker='o',markersize='6.0')
    plt.ylabel("|atomic forces| (eV/$\AA$)")
    plt.xlabel(" Distance of the host atoms from the muon ($\AA$)")
    return


def SCconv_withscf_forces(pw_out,mu_atomic_num):   #takes in  ASE Atom structure data format with forces and Integer atomic number of the muon
    #pw_out is ASE Atom structure data format with forces
    #mu_atomic_num is Integer atomic number of the muon, !!!MAYBE USING ATOMIC SYMBOL IS BETTER.
    conv_thresh=1e-03*Rydberg/Bohr  #from au to eV/A
    # 
    #
    mu_idd=[atom.index for atom in pw_out if atom.number == mu_atomic_num]   #Get muon index
    if len(mu_idd)>1:
        raise SystemExit("!!!Check provided muon atomic number has more than one muon in the structure. Exiting....") from None
    mu_id=mu_idd[0]
    mu_symb=pw_out.symbols[mu_id]    #mu symbol
    #
    #
    atm_forces_vec=pw_out.get_forces()  # atomic forces vector form in eV/A
    atm_forces_mag=[np.sqrt(x.dot(x)) for x in atm_forces_vec]  #magnitude of atomic forces
    #
    ##mu_force=atm_forces_mag.pop(mu_id) #remove mu force from the list, we don't want to mess with the index of stored forces
    no_mu_atm_forces_mag1=[v for i,v in enumerate(atm_forces_mag) if i !=mu_id] # forces without that of muon
    if min(no_mu_atm_forces_mag1) < conv_thresh: 
        print("First SC convergence Criteria achieved", " min force in au",min(no_mu_atm_forces_mag1)*Bohr/Rydberg)
        cond = True
    else:
        print("!!!First SC convergence Criteria  NOT achieved", " min force in au",min(no_mu_atm_forces_mag1)*Bohr/Rydberg)
        cond=False
    #
    #
    #atm_indxes=[atom.index for atom in pw_out if atom.index !=  mu_id]  #atom indeces except that of muon
    atm_indxes=[atom.index for atom in pw_out]  #atom indeces
    atm_num=len(pw_out.numbers)                 # number of atoms
    atm_dist=pw_out.get_distances( mu_id, atm_indxes, mic=True, vector=False) # minimum image distance from muon in Ang
    #
    #
    specie_num=len(set(pw_out.numbers))  # number of species
    specie_set=set(pw_out.get_chemical_symbols()) #set of specie
    specie_set = [i for i in specie_set if i !=mu_symb] #remove mu symbol from the iteration list
    cond2=[]
    for i in range(0,specie_num-1):
        specie_index=[atom.index for atom in pw_out if atom.symbol == specie_set[i]]
        if len(specie_index)> 5:
            #
            #specie_dist=[atm_dist[i] for i in specie_index]
            #specie_force=[atm_forces_mag[i] for i in specie_index]
            #
            # from trend, seems forces above 0.06 au are the ones closest to the muon, or its better to compare atm_dist and atm_forces_mag
            specie_dist=[atm_dist[i] for i in specie_index if atm_forces_mag[i] < 0.06*Rydberg/Bohr] 
            specie_force=[atm_forces_mag[i] for i in specie_index if atm_forces_mag[i] < 0.06*Rydberg/Bohr]
            par,cov=fit_curve(exp_fnc,specie_dist, specie_force)
            #
            #stder= np.sqrt(np.diag(cov))   
            #if stder[0] > par[0] or stder[1] > par[1]: #fit  and data check- worse case scenario!!! better conditions?
            #    raise SystemExit("!!!Check data and fit,Maybe the data does not decay exponentialy. Exiting....") from None
            #    
            min_conv_dist=min_SCconv_dist(conv_thresh,par[0],par[1]) #func to find min distance for conv
            #
            # plotting
            color=['purple','orange','green','red','black','blue','yellow']
            plot_forces_vs_dist_frm_mu(specie_dist,specie_force,par,specie_set[i],color[i])
            #
            #
            print("Max mu-atom distance in the SC is",max(specie_dist),", Min distance required for conv is",min_conv_dist)
            if max(specie_dist)>= min_conv_dist:
                print('For specie',specie_set[i],", Second SC convergence Criteria achieved")
                #cond2 = True
                cond2.append(True)
            else:
                print('!!!For specie',specie_set[i],", Second SC convergence NOT achieved","SC with min",min_conv_dist,"Ang required")
                cond2.append(False)
        else:
            print('!!!For specie',specie_set[i],', the current SC size is NOT sufficient for convergence checks')
            cond2.append(False)
    
    return cond, cond2

if __name__ == "__main__":
    #pw_out=read("cif_pwo_files/LiF_scf_1x1x1.out",index=-1,format='espresso-out')
    pw_out=read("cif_pwo_files/LiF_scf_2x2x2.out",index=-1,format='espresso-out')
    #pw_out=read("cif_pwo_files/LiF_scf_3x3x3.out",index=-1,format='espresso-out')
    #pw_out=read("cif_pwo_files/mnbite2_relax_1+.out",index=0,format='espresso-out')
    #pw_out=read("cif_pwo_files/fe.relax2x2Hoct.out",index=0,format='espresso-out')
    #pw_out=read("cif_pwo_files/fe.relax2x2Htet.out",index=0,format='espresso-out')
    #pw_out=read("cif_pwo_files/fe.relax3x3Hoct.out",index=0,format='espresso-out')
    #pw_out=read("cif_pwo_files/fe.relax3x3Htet.out",index=0,format='espresso-out')
    cond,cond2=SCconv_withscf_forces(pw_out,1) 
    plt.show()