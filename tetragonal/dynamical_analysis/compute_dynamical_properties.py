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
params["nstates"]        = 76 # total number of electronic states
params["init_times"]     = [600]
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
params["nstates"]        = 118 # total number of electronic states
params["init_times"]     = [600]
params["active_space"]   = list(range(0,76))#params["nstates"])) # indexing is from 0!
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
    returns: a list of energies vs time for all states in hvib
    """
    nsteps   = len(hvib)
    nstates  = hvib[0].num_of_rows
    energies = []
    for state in range( nstates ):
        energies.append( [] )
        for step in range( nsteps ):
            energies[ state ].append( hvib[ step ].get( state, state ).real - hvib[ step ].get( 0, 0 ).real )
    return np.array( energies )


def compute_tnacs( hvib ):
    """
    Computes the time-averaged nonadiabatic couplings for a given hvib.
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
    returns: a matrix of time-averaged nonadiabatic couplings between all electronic states
             in meV
    """
    nstates  = hvib[0].num_of_rows
    #nacs = data_stat.cmat_stat2( hvib, 2)
    nacs = data_stat.cmat_stat2( hvib, 1)
    mb_tnacs = []
    for i in range( nstates ):
        mb_tnacs.append( [] )
        for j in range( nstates ):
            mb_tnacs[i].append( nacs.get(i,j).imag * 1000.0 / units.ev2Ha )
    return np.array( mb_tnacs )



#####################
# 3. Divide up the data (obtained from 1) into sub-trajectories. These are to be consdiered our independent nuclear sub-trajectories.
#    Here though, we just consider 1 subtrajectory, which is the entire trajectory. One could however consider more
#    The first zero below corresponds to the actual MD trajectory index (from zero). We have only 1 MD trajectory computed,
#    so this number is zero. The second number corresponds to the start time of the MD trajectory. Shown (but commented) 
#    is how one would consdier a "sub" trajectory that is the Nth MD trajectory starting from time t
#    Each trajectory, starting at time "t" is considered for length "subtraj_len"  

subtraj_time_info = [
                      [0, 0], #[N, t],
                    ]

nsubtrajs   = len(subtraj_time_info)
subtraj_len = params["nfiles"] # steps
params["dt"] = 1.0 * units.fs2au

hvib_mb_subtrajs = []
hvib_mixed_sd_subtrajs = []
hvib_elec_sd_subtrajs = []
hvib_hole_sd_subtrajs = []

nstates_mb       = hvib_mb[0][0].num_of_rows
nstates_mixed_sd = hvib_mixed_sd[0][0].num_of_rows

for i in range( nsubtrajs ):
    hvib_mb_subtrajs.append( [] )
    hvib_mixed_sd_subtrajs.append( [] )


# This is a function that is python mp.pool parallelized.

def myfunc( subtraj ):

    nuclear_traj = subtraj_time_info[subtraj][0]
    start_time   = subtraj_time_info[subtraj][1]

    print("\nsubtraj", subtraj)
    print("nuclear trajectory", nuclear_traj)
    print("start time =", start_time)
    print("end time   =", start_time + subtraj_len)

    md_time = [ i for i in range( subtraj_len ) ]
    md_time = np.array( md_time ) * params["dt"] * units.au2fs

    # Obtain the subtrajectory
    for i in range(start_time, start_time + subtraj_len):
        hvib_mb_subtrajs[ subtraj ].append( hvib_mb[nuclear_traj][i] )
        hvib_mixed_sd_subtrajs[ subtraj ].append( hvib_mixed_sd[nuclear_traj][i] )

    # Compute properties
    mb_energies = compute_state_energies_vs_time( hvib_mb_subtrajs[ subtraj ] )
    mb_tnacs    = compute_tnacs( hvib_mb_subtrajs[ subtraj ] )
    mb_tau, mb_rates = decoherence_times.decoherence_times( hvib_mb_subtrajs[ subtraj ] )

    mixed_sd_energies = compute_state_energies_vs_time( hvib_mixed_sd_subtrajs[ subtraj ] )
    mixed_sd_tnacs    = compute_tnacs( hvib_mixed_sd_subtrajs[ subtraj ] )
    mixed_sd_tau, mixed_sd_rates = decoherence_times.decoherence_times( hvib_mixed_sd_subtrajs[ subtraj ] )

    print("\nTesting np.mean np.std for mb case")
    for state in range( nstates_mb ):
        s1_energy = np.mean(mb_energies[1])
        print("State", state, np.mean(mb_energies[state]-s1_energy)*units.au2ev)

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    #plt.title('Tetragonal, MB', fontsize=10) # traj'+str(subtraj)+'', fontsize=10)
    plt.title('Tetragonal', fontsize=10) # traj'+str(subtraj)+'', fontsize=10)
    plt.xlabel('Time, fs',   fontsize=10)
    plt.ylabel('Energy, eV', fontsize=10)
    plt.ylim(1.5,2.3)
    plt.yticks([2,2.3])
    #plt.xticks([0,10,20,30,40,50])
    for mb_index in range( nstates_mb ):
        plt.plot(md_time, mb_energies[mb_index]*units.au2ev, label="", linewidth=1, color="black")
    plt.tight_layout()
    #plt.legend(fontsize=8, ncol=3, loc="lower left")
    plt.savefig("CsPbI3_222tetra_mb_Energies_"+str(subtraj)+".png", dpi=300)

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title('Tetragonal, MB', fontsize=10) # traj'+str(subtraj)+'', fontsize=10)
    #plt.xlabel("Many-body basis")
    #plt.ylabel("Many-body basis")
    plt.xlabel("State Index")
    plt.ylabel("State Index")
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlim(0,nstates_mb-1)
    plt.ylim(0,nstates_mb-1)
    plt.xticks(range(0,nstates_mb,25))
    plt.yticks(range(0,nstates_mb,25))
    plt.imshow( mb_tnacs, cmap='hot', interpolation='nearest', vmin=0, vmax=30)
    cb = plt.colorbar(label="meV")
    cb.ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig("CsPbI3_222tetra_mb_tnacs.png")

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title('Tetragonal, SP', fontsize=10) # traj'+str(subtraj)+'', fontsize=10)
    plt.xlabel('Time, fs',   fontsize=10)
    plt.ylabel('Energy, eV', fontsize=10)
    plt.ylim(1,3.6)
    plt.yticks([1,2,3])
    for sd_index in range( nstates_mixed_sd ):
        plt.plot(md_time, mixed_sd_energies[sd_index]*units.au2ev, label="", linewidth=1)
    plt.tight_layout()
    #plt.legend(fontsize=8, ncol=3, loc="lower left")
    plt.savefig("CsPbI3_222tetra_mixed_sd_Energies_"+str(subtraj)+".png", dpi=300)

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title('Tetragonal, SP', fontsize=10) # traj'+str(subtraj)+'', fontsize=10)
    plt.xlabel("Slater determinant basis")
    plt.ylabel("Slater determinant basis")
    plt.xlabel("State Index")
    plt.ylabel("State Index")
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlim(0,nstates_mixed_sd-1)
    plt.ylim(0,nstates_mixed_sd-1)
    plt.xticks(range(0,nstates_mixed_sd,25))
    plt.yticks(range(0,nstates_mixed_sd,25))
    plt.imshow( mixed_sd_tnacs, cmap='hot', interpolation='nearest', vmin=0, vmax=30)
    cb = plt.colorbar(label="meV")
    cb.ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig("CsPbI3_222tetra_mixed_sd_tnacs.png")




pool = mp.Pool(nsubtrajs)
pool.map( myfunc, list(range(nsubtrajs)) )
pool.close()
pool.join()


