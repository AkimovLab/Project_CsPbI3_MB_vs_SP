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
import matplotlib.pyplot as plt
import os

colors = {}
colors.update({"1": '#000000'})  # Black 
colors.update({"2": '#000099'})  # Blue  
colors.update({"3": '#006400'})  # Green 
colors.update({"4": '#990000'})  # Red   
colors.update({"5": '#8B008B'})  # Purple
colors.update({"6": '#FF8C00'})  # Orange
colors.update({"9": '#4d4d4d'})  # Gray  
color_index = ["1","2","3","4","5","6","9"]





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







#####################
# 2. Divide up the data (obtained from 1) into sub-trajectories. These are to be consdiered our independent nuclear sub-trajectories.
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

for subtraj in range( nsubtrajs ):

    hvib_mb_subtrajs.append( [] )
    hvib_mixed_sd_subtrajs.append( [] )

    nuclear_traj = subtraj_time_info[subtraj][0]
    start_time   = subtraj_time_info[subtraj][1]

    print("nuclear trajectory", nuclear_traj)
    print("start time =", start_time)
    print("end time   =", start_time + subtraj_len)

    md_time = [ i for i in range( subtraj_len ) ]
    md_time = np.array( md_time ) * params["dt"] * units.au2fs

    # Obtain the subtrajectory
    for i in range(start_time, start_time + subtraj_len):
        hvib_mb_subtrajs[ subtraj ].append( hvib_mb[nuclear_traj][i] )
        hvib_mixed_sd_subtrajs[ subtraj ].append( hvib_mixed_sd[nuclear_traj][i] )

    ###### Obtain all nacs for all pairs states - put into one giant list
    # define some lists and variables
    nac_mb, nac_sd = [], []
    states_mb = range( hvib_mb[0][0].num_of_rows )
    states_sd = range( hvib_mixed_sd[0][0].num_of_rows )
    nac_cutoff = 0   # meV for zero, this means we consider all possible values for nacs

    for t in range(start_time, start_time + subtraj_len):

        # MB
        for i in states_mb:
            for j in states_mb:    
                if j != i:
                    nac_mb_ij = abs( hvib_mb[0][t].get(i,j).imag ) * 1000.0 * units.au2ev
                    if nac_mb_ij > nac_cutoff:
                        x_mb = MATRIX(1,1)
                        x_mb.set(0, 0, nac_mb_ij )
                        nac_mb.append( x_mb )

        # SD
        for i in states_sd:
            for j in states_sd:
                if j != i:
                    nac_sd_ij = abs( hvib_mixed_sd[0][t].get(i,j).imag ) * 1000.0 * units.au2ev
                    if nac_sd_ij > nac_cutoff:
                        x_sd = MATRIX(1,1)
                        x_sd.set(0, 0, nac_sd_ij )
                        nac_sd.append( x_sd )

    #==========================
    # Sanity check
    nac_mb[0].show_matrix()
    print( hvib_mb[0][0].get(0,1) * 1000.0 * units.au2ev )
    print("len(nac_mb) is :",len(nac_mb))
    print("len(nac_sd) is :",len(nac_sd))
    # End of Sanity check
    #==========================

    # For MB
    bin_supp_mb, dens_mb, cum_mb = data_stat.cmat_distrib( nac_mb, 0, 0, 0, 0, 50, 0.1)
    # For SD
    bin_supp_sd, dens_sd, cum_sd = data_stat.cmat_distrib( nac_sd, 0, 0, 0, 0, 50, 0.1)




fs = 12
plt.figure(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
plt.subplot(1,1,1)
plt.title('Tetragonal', fontsize=fs)
plt.xlabel('|NAC|, meV', fontsize=fs)
plt.ylabel('PD, 1/meV',  fontsize=fs)
plt.yticks(fontsize=fs)
plt.xticks(fontsize=fs)
plt.xlim(0,5)
plt.plot( bin_supp_mb, dens_mb, label="MB", linewidth=2.0, color = "blue")
plt.plot( bin_supp_sd, dens_sd, label="SP", linewidth=2.0, color = "red")
plt.legend(fontsize=12, loc='upper right')
plt.tight_layout()
plt.savefig("Tetragonal_MB_vs_SD_1.png", dpi=600)

plt.figure(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
plt.subplot(1,1,1)
plt.title('Tetragonal', fontsize=fs)
plt.xlabel('|NAC|, meV', fontsize=fs)
plt.ylabel('PD, 1/meV',  fontsize=fs)
plt.xlim(5,20)
plt.ylim(0,0.02)
plt.plot( bin_supp_mb, dens_mb, label="MB", linewidth=2.0, color = "blue")
plt.plot( bin_supp_sd, dens_sd, label="SP", linewidth=2.0, color = "red")
plt.yticks(ticks=[0, 0.005, 0.01, 0.015, 0.02],labels=['0','5','10','15',20],fontsize=fs)
plt.text(5.5, 0.021, 'x10$^{-2}$', fontsize=12)
plt.legend('',frameon=False)
plt.tight_layout()
plt.savefig("Tetragonal_MB_vs_SD_2.png", dpi=600)


plt.figure(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
plt.subplot(1,1,1)
plt.title('Tetragonal', fontsize=fs)
plt.xlabel('|NAC|, meV', fontsize=fs)
plt.ylabel('PD, 1/meV',  fontsize=fs)
plt.xlim(20,50)
plt.ylim(0,0.002)
plt.plot( bin_supp_mb, dens_mb, label="MB", linewidth=2.0, color = "blue")
plt.plot( bin_supp_sd, dens_sd, label="SP", linewidth=2.0, color = "red")
plt.yticks(ticks=[0, 0.0005, 0.001, 0.0015, 0.0020],labels=['0','5','10','15',20],fontsize=fs)
plt.text(21, 0.0021, 'x10$^{-3}$', fontsize=12)
plt.legend('',frameon=False)
plt.tight_layout()
plt.savefig("Tetragonal_MB_vs_SD_3.png", dpi=600)

