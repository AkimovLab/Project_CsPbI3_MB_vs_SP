import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize
from libra_py import data_conv
from libra_py import data_stat
import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import os





####################
# 1. Read / Get the Hvib files from step3/res_mb_sp
# For the CI (MB) case 
print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../step3/res_mb_sp/")
params["Hvib_re_prefix"] = "Hvib_ci_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_ci_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 1800
params["nstates"]        = 151 # total number of electronic states
params["init_times"]     = [0]
params["active_space"]   = list(range(params["nstates"])) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_mb = step4.get_Hvib2(params)
hvib_mb[0][-1].show_matrix()
print ("Length of hvib_mb is: ", len(hvib_mb[0]))
#sys.exit(0)

# For the SD (SP) case
print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../step3/res_mb_sp/")
params["Hvib_re_prefix"] = "Hvib_sd_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_sd_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 1800
params["nstates"]        = 229 # total number of electronic states
params["init_times"]     = [0]
params["active_space"]   = list(range(0,151))#params["nstates"])) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_mixed_sd = step4.get_Hvib2(params)
hvib_mixed_sd[0][-1].show_matrix()
print ("Length of hvib_mixed_sd is: ", len(hvib_mixed_sd[0]))
#sys.exit(0)





####################
# 2. Define some helper functions. Read each function for more detail

def compute_state_energies_vs_time( hvib ):
    """
    Computes the states energies vs time for a given hvib.
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
    returns: a list of CMATRIX objects that hold the energies vs time for all states in hvib
    """
    nsteps   = len(hvib)
    nstates  = hvib[0].num_of_rows
    energies = []

    print("nsteps", nsteps)
    print("nstates", nstates)
    for step in range( nsteps ):
        # For every step, subtract the ground state energy
        gs_en = CMATRIX( nstates, nstates )
        for state in range( nstates ):
            gs_en.set( state, state, hvib[step].get(0,0) )
        energies.append( ( hvib[step] - gs_en ) * units.au2ev )
    return energies








# Some needed initialization
dens_mb_ii = []
dens_mixed_sd_ii = []
nstates_mb = hvib_mb[0][0].num_of_rows
nstates_mixed_sd = hvib_mixed_sd[0][0].num_of_rows
print("nstates", nstates_mb, nstates_mixed_sd)

# Compute the state energies vs. time. These are lists of CMATRIX objects
mb_energies = compute_state_energies_vs_time( hvib_mb[0] )
mixed_sd_energies = compute_state_energies_vs_time( hvib_mixed_sd[0] )

# Compute the probability density distributions for each state, then take their average
deV = 0.01
for i in range(1,nstates_mb):
    bin_supp, dens, cum = data_stat.cmat_distrib( mb_energies, i, i, 0, 0.0, 4, deV)    
    dens_mb_ii.append( np.array(dens)*deV )
mb_dens = sum(dens_mb_ii)

for i in range(1,nstates_mixed_sd):
    bin_supp, dens, cum = data_stat.cmat_distrib( mixed_sd_energies, i, i, 0, 0.0, 4, deV)
    dens_mixed_sd_ii.append( np.array(dens)*deV )
mixed_sd_dens = sum(dens_mixed_sd_ii)


# Ploting below
plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
plt.subplot(1,1,1)
plt.title('Cubic', fontsize=10) # traj'+str(subtraj)+'', fontsize=10)
plt.ylabel('DOS, 1/eV',  fontsize=10)
plt.xlabel('Energy, eV', fontsize=10)
plt.xlim(1,4)
plt.xticks([1,2,3,4])
#plt.yticks([])
plt.plot(bin_supp, mixed_sd_dens, color="red", label="SD")
plt.plot(bin_supp, mb_dens, color="blue", label="MB")
plt.legend()
plt.tight_layout()
plt.savefig("DOS.png")
sys.exit(0)


