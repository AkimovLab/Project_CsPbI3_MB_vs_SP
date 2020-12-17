import os
import sys
import time
import numpy as np
import multiprocessing as mp
import glob
import matplotlib.pyplot as plt

from libra_py import pdos

####################################################################################################################################################
#============================= Main part starts from here and we use the function above combined with multiprocessing =============================#

#pdos_type = "orbital_resolved"
pdos_type = "atom_resolved"
#thermal=True
thermal=False

if pdos_type == "orbital_resolved":

    #===================== Orbital resolved columns
    #                               Cs, total            Cs, s             Cs, p             Cs, d
    angular_momentum_cols = [ [ list(range(3,12)), list(range(3,4)), list(range(4,7)), list(range(7,12)) ],
    #                               Pb, total            Pb, s             Pb, p             Pb, d
                              [ list(range(3,12)), list(range(3,4)), list(range(4,7)), list(range(7,12)) ],
    #                               I, total             I, s              I, p              I, d
                              [ list(range(3,12)), list(range(3,4)), list(range(4,7)), list(range(7,12)) ]
                            ] 

    # This is for orbital resolved
    labels = ['Cs, s','Cs, p','Cs, d', 'Pb, s','Pb, p','Pb, d', 'I, s','I, p','I, d']
    # Colors, the color orders are based on the labels, Here are chosen for elements
    colors = ['red','blue','green','orange','purple','pink','gold','cyan','brown']
    outname = 'orbital'

elif pdos_type == "atom_resolved":

    #===================== Atom resolved columns
    #                           Cs, total             Cs, total
    angular_momentum_cols = [ [ list(range(3,12)), list(range(3,12)) ],
    #                           Pb, total             Pb, total
                              [ list(range(3,12)), list(range(3,12)) ],  
    #                           I, total              I, total
                              [ list(range(3,12)), list(range(3,12)) ] 
                            ]

    # This is for atom resolved
    labels = ['Cs','Pb', 'I']
    # This row is for atom resolved
    colors = ['blue','green', 'purple']
    outname = 'atom'





#===================== Other inputs
# The time step is MD, For static calculations we only use 0
time_step = 0
sigma = 0.1
coef = 1
npoints = 2000 # number of points for the grid mesh vector
energy_conversion = 27.211386 # Hartree to eV
nprocs = 24 # number of processors

t1 = time.time()
# create the pool of processors
pool = mp.Pool(nprocs)

# This variable will contain the convolved PDOS for each element and each angular momentum list
Total_results = []
homos = []
#========================== Convolution
# Here we change the order,
# In fact the k1 is Cs, k2 is Pb, k3 is I
for k in [1,2,3]:

    # Create an empty list for summation of the convolved PDOS
    # This is used for Total DOS
    dos_summation_angular = []
    energy_grid_ave = []
    
    # Define a zero vector to sum the convolved PDOS for each angular momentum list
    # It is based on the number of points
    zero_vec = np.zeros((npoints))
    # Now for each angular momentum append the zero vector
    # The same for enegy grid average since we want the energy grid average
    for i in range(len(angular_momentum_cols[k-1])):
        dos_summation_angular.append(zero_vec)
        energy_grid_ave.append(zero_vec)

    # Create the variables for pool.starmap
    vars_for_pool = []
    # Find all the pdos files of an element in all_pdosfiles for the first trajectory

    if thermal == True:
        DOS_files1 = glob.glob('all_pdosfiles/*k%d-1.pdos'%k)
    else:
        DOS_files1 = glob.glob('../tddft/*k%d-1.pdos'%k)

    for DOS_file in DOS_files1:
        params = {}
        params["cp2k_pdos_file"] = DOS_file
        params["time_step"] = time_step
        params["sigma"] = sigma
        params["coef"] = coef
        params["npoints"] = npoints
        params["energy_conversion"] = energy_conversion
        params["angular_momentum_cols"] = list(angular_momentum_cols[k-1])
        vars_for_pool.append(params)

    results_for_k = pool.map(pdos.convolve_cp2k_pdos, vars_for_pool)

    # We initialize all the homos average: homos_ave
    # This variable is the same for each element so it will be repeated but doesn't 
    # change the results
    homos_ave = 0
    # Now we need to take the average of them (The same is for average HOMO eergy level
    # and the energy grid as well)
    # Energy grid is the 0th element, convolved DOS is the 1st element, HOMO energy levels is the 2nd one
    # 1. We add the convolved ones to dos_summation_angular
    # Here we start to take the averages by summing the results.
    for i in range(len(DOS_files1)):
        energy_grid_ave       += results_for_k[i][0] # First element is the energy grid output by the function
        dos_summation_angular += results_for_k[i][1] # Second element is the convolved DOS
        homos_ave             += results_for_k[i][2] # Third element is the HOMO energy level
        
        
    # 2. We take the average for that by dividing by len(DOS_files)
    # 3. We append it to Total_results
    Total_results.append(dos_summation_angular/len(DOS_files1))
    # Uncomment only if you use two trajectories
    #Total_results.append(dos_summation_angular/(len(DOS_files1)+len(DOS_files2)))
    energy_grid_ave /= len(DOS_files1)
    # Uncomment only if you use two trajectories
    #energy_grid_ave /= (len(DOS_files1)+len(DOS_files2))
    homos_ave /= len(DOS_files1)
    # Uncomment only if you use two trajectories
    #homos_ave /= (len(DOS_files1)+len(DOS_files2))
    

# Close the pool
pool.close()
pool.join()
#================================== End of convolution



#================================== Total density
# We first compute the total density of states through the first computed angular momentum 
# for [3,12] as is shown above - So we only choose the 0 index
# Make a zero vector of npoints
total_density = np.zeros((npoints))
# i here is each element
for i in range(len(Total_results)):
    # Sum the total by Total_results of zero index angular momentum column for each element
    total_density += Total_results[i][0]

#============================== Plotting 
figure = plt.figure(num=None, figsize=(3.21, 2.41), dpi=1200, edgecolor='black', frameon=True)


# Plot the total density by black color
plt.plot(energy_grid_ave[0]-homos_ave,total_density,label='Total',color='black', linewidth=2.0)


# set up a counter for labels
c = 0
for i in range(len(Total_results)):
    for j in range(1,len(Total_results[i])):
        plt.plot(energy_grid_ave[0]-homos_ave,Total_results[i][j],label=labels[c],color=colors[c], linewidth=2)
        c += 1



plt.xlim(-4,6)
plt.ylim(0,20)
plt.legend(fontsize=6.75, ncol=1, loc='upper right')
plt.xlabel('Energy, eV',fontsize=12)
plt.ylabel('DOS, 1/eV',fontsize=12)

if thermal == True:
    plt.title("Tetragonal, 300 K",fontsize=12)
else:
    plt.title("Tetragonal, 0 K",fontsize=12)

plt.tight_layout()
plt.show()

if thermal == True:
    if pdos_type == "orbital_resolved":
        plt.savefig('orbital_DOS_300K_'+outname+'_average.png', dpi=300)
    elif pdos_type == "atom_resolved":
        plt.savefig('atom_DOS_300K_'+outname+'_average.png', dpi=300)

elif thermal == False:
    if pdos_type == "orbital_resolved":
        plt.savefig('orbital_DOS_0K_'+outname+'.png', dpi=300)
    elif pdos_type == "atom_resolved":
        plt.savefig('atom_DOS_0K_'+outname+'.png', dpi=300)

print("Total Elapsed time:",time.time()-t1)

