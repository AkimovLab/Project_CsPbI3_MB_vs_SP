import os
import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize
from libra_py import data_conv
from libra_py import data_read
from libra_py import data_stat
import libra_py.workflows.nbra.mapping as mapping
import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.step3 as step3
import libra_py.workflows.nbra.step2_many_body as step2_many_body
import libra_py.workflows.nbra.step2_analysis as step2_analysis
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

#######################################################################################
# Helper functions
def compute_overlaps_in_parallel( step ):
    """
    Function used for the making the computation of overlaps parallel via python multiprocessing 
    """
    s_sd  = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step], S_ks[0][step],  use_minimal=False )
    s_sd  = data_conv.MATRIX2nparray(s_sd)
    return s_sd

def compute_toverlaps_in_parallel( step ):
    """
    Function used for the making the computation of time-overlaps parallel via python multiprocessing 
    """
    st_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step+1], St_ks[0][step], use_minimal=False )
    st_sd = data_conv.MATRIX2nparray(st_sd)
    return st_sd



#######################################################################################
"""
 1. Read the files that have the energies, overlaps, and time-overlap matricies in the 
    Kohn-Sham basis: E_ks, S_ks, and St_ks.
"""
#######################################################################################

os.system("rm -rf res_mb_sp")
res_dir = "res_mb_sp"
os.mkdir(res_dir)

path = os.getcwd()
res_dir_mb   = path+"/../step2/res/"
data_dim     = 132 # rows in E_ks
active_space = range(0,int(data_dim/2)) # alpha channel only here #range(data_dim)
istep        = 600  # initial step
fstep        = 2401 # final step + 1  
dt           = 1.0*units.fs2au

params = { "data_set_paths" : [res_dir_mb], "data_dim":data_dim, "active_space":active_space, "isnap":istep,  "fsnap":fstep }
# Fetching E_ks
params.update({ "data_re_prefix" : "E_ks_",  "data_re_suffix" : "_re", "data_im_prefix" : "E_ks_",  "data_im_suffix" : "_im"  } )
E_ks = data_read.get_data_sets(params)
E_ks[0][-1].show_matrix()

# Fetching S_ks
params.update({ "data_re_prefix" : "S_ks_", "data_re_suffix" : "_re", "data_im_prefix" : "S_ks_", "data_im_suffix" : "_im"  } )
S_ks = data_read.get_data_sets(params)

# Fetching St_ks
params.update({ "data_re_prefix" : "St_ks_", "data_re_suffix" : "_re", "data_im_prefix" : "St_ks_", "data_im_suffix" : "_im"  } )
St_ks = data_read.get_data_sets(params)
St_ks[0][-2].show_matrix()
#sys.exit(0)

step3.apply_orthonormalization_general( S_ks[0], St_ks[0] )
#step3.apply_phase_correction( St_ks[0] )
S_ks[0][-1].show_matrix()
#sys.exit(0)







#######################################################################################
"""
 2. Read the TDDFT output files and get all the needed information 
    2.1. Set parameters
    2.2. Get the information from the TDDFT calculations
    2.3. Reindex the single-particle excitations into a format expected by Libra
    2.4. Order / sort the single-particle excitations at each timestep by energy or identity
"""
#######################################################################################

# 2.1. Set parameters
num_excited_states = 75
min_band_ks   = 515
max_band_ks   = 580
ks_homo_index = 544
ks_orbital_indicies = list( range( min_band_ks, max_band_ks + 1 ) )
params["ks_orbital_indicies"] = ks_orbital_indicies
params["logfile_directory"]   = "../excitation_analysis/all_logfiles"
params["es_software"]         = "cp2k"
params["isUKS"]               = 0
params["number_of_states"]    = num_excited_states 
params["tolerance"]           = 0.0
params["start_time"]  = istep
params["finish_time"] = fstep

# 2.2. Get the information from the TDDFT calculations
S_sd, St_sd, sd_basis_states_unique, ci_basis_states, ci_coefficients, ci_energies, spin_components = step2_analysis.get_step2_mb_sp_properties( params )

# 2.3. Reindex the single-particle excitations into a format expected by Libra
sd_states_reindexed = step2_many_body.reindex_cp2k_sd_states( ks_homo_index, ks_orbital_indicies, sd_basis_states_unique, sd_format=2 )

# 2.4. Order / sort the single-particle excitations at each timestep by energy or identity
E_sd, sd_states_unique_sorted, sd_states_reindexed_sorted, reindex_nsteps = step3.sort_unique_SD_basis( E_ks[0], sd_basis_states_unique, sd_states_reindexed, istep, fstep, sorting_type="energy" )
#sys.exit(0)







#######################################################################################
"""
 3. Compute overlaps and time-overlaps between the single-particle excitations 
    3.1. Compute overlaps in the basis of single-particle excitations
    3.2. Compute time-overlaps in the basis of single-particle excitations
    3.3. Convert the data types to Libra's CMATRIX
"""
#######################################################################################

# 3.1. Overlaps
pool = mp.Pool( 24 )
tmp_s_sd = pool.map( compute_overlaps_in_parallel, list( range( fstep - istep ) ) )
pool.close()
pool.join()

# 3.2. Time-overlaps
pool = mp.Pool( 24 )
tmp_st_sd = pool.map( compute_toverlaps_in_parallel, list( range( fstep - istep - 1 ) ) )
pool.close()
pool.join()

# 3.3. Convert the data types to Libra's CMATRIX
for step in range( fstep - istep ):
    S_sd.append(  data_conv.nparray2CMATRIX( tmp_s_sd[step]  ) )
for step in range( fstep - istep - 1 ):
    St_sd.append( data_conv.nparray2CMATRIX( tmp_st_sd[step] ) )
#sys.exit(0)






#######################################################################################
"""
 4. 
    4.1. Take the list of excitation energies at each timestep, and compute the midpoints
    4.2. Make the transformation matrix from the single-particle to many-body basis at each 
         time step. 
    4.3. Transform from single-particle to many-body (CI) basis
    4.4. Apply orthonormalization to the many-body basis, apply state reordering, 
    and apply phase corrections
    4.5. Output many-body basis overlaps and time-overlaps to the res directory
    4.6. Make the Hvib in the many-body basis
"""
#######################################################################################


# 4.1. Take the list of excitation energies at each timestep, and compute the midpoints 
ci_midpoint_energies = step3.compute_ci_energies_midpoint( ci_energies, num_excited_states, istep, fstep )
#sys.exit(0)

# 4.2. Make the transformation matrix from the single-particle to many-body basis at each time step
SD2CI = step3.make_T_matricies( ci_coefficients, ci_basis_states, spin_components, sd_states_unique_sorted, num_excited_states, istep, fstep, res_dir )
#sys.exit(0)

# 4.3. Transform from single-particle to many-body (CI) basis
S_ci, St_ci  = [], []
for step in range( fstep - istep ):
    s_ci  = SD2CI[step].H() * S_sd[step]  * SD2CI[step]
    S_ci.append(  s_ci  )
for step in range( fstep - istep - 1 ):
    st_ci = SD2CI[step].H() * St_sd[step] * SD2CI[step+1]
    St_ci.append( st_ci )

# 4.4. Apply orthonormalization to the many-body basis, apply state reordering, apply phase corrections
step3.apply_orthonormalization_general( S_ci, St_ci )
params2 = {"do_state_reordering":2, "state_reordering_alpha":0}
step3.apply_state_reordering_general( St_ci, ci_midpoint_energies, params2 )
step3.apply_phase_correction_general( St_ci )

# 4.5. Output many-body basis overlaps and time-overlaps to the res directory
print("Outputting the CI data to the res directory..." )
for step in range( fstep - istep ):
    S_ci[step].real().show_matrix("%s/S_ci_%d_re" % (res_dir, int(istep+step)))
for step in range( fstep - istep - 1 ):
    St_ci[step].real().show_matrix("%s/St_ci_%d_re" % (res_dir, int(istep+step)))

# 4.6. Make the Hvib in the many-body basis
for step in range( fstep - istep - 1 ):
    ci_nacs = (  0.5j / dt ) * CMATRIX ( ( St_ci[step] - St_ci[step].H() ).real() )    
    ci_hvib = ci_midpoint_energies[step] - ci_nacs
    ci_hvib.real().show_matrix("%s/Hvib_ci_%d_re" % (res_dir, int( istep+step )))
    ci_hvib.imag().show_matrix("%s/Hvib_ci_%d_im" % (res_dir, int( istep+step )))








#######################################################################################
"""
 5. Make and output the Hvib in the basis of single-particle excitations
    5.1. Compute single-particle excitation energies at the midpoints
    5.2. Apply orthonormalization to the many-body basis, apply state reordering, 
    and apply phase corrections
    5.3. Make the hvib in the basis of single-particle excitations
"""
#######################################################################################

# 5.1. Compute single-particle excitation energies at the midpoints
sd_midpoint_energies = []
for step in range( fstep - istep - 1 ):
    sd_midpoint_energy = 0.5 * ( E_sd[step] + E_sd[step+1] )
    sd_midpoint_energies.append( sd_midpoint_energy )
print("\nNormalize SD basis before output")

# 5.2. Apply orthonormalization to the many-body basis, apply state reordering, and apply phase corrections
step3.apply_orthonormalization_general( S_sd, St_sd )
params2 = {"do_state_reordering":2, "state_reordering_alpha":0}
step3.apply_state_reordering_general( St_sd, sd_midpoint_energies, params2 )
step3.apply_phase_correction_general( St_sd )

# 5.3. Make the Hvib in the basis of single-particle excitations
for step in range( fstep - istep - 1 ):
    sd_nacs = (  0.5j / dt ) * CMATRIX ( ( St_sd[step] - St_sd[step].H() ).real() )
    sd_hvib = sd_midpoint_energies[step] - sd_nacs
    sd_hvib.real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir, int( istep+step )))
    sd_hvib.imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir, int( istep+step )))
# END


