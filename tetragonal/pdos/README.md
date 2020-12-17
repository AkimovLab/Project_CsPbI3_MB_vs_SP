## pDOS

This file contains the scripts needed to compute the projected density of states calculations using CP2K and Libra 

In compute_pdos.py, set "thermal" = False when computing the projected density of states for the optimized system. Set set "thermal" = False when computing the projected density of states for the thermally averaged system.

For the thermally averaged system, you will need a folder called "all_pdosfiles" present in this directy, and it should contain the output (.pdos) file from each of the KS-DFT calculations from cp2k, you can generate this directory by running the following commands (assuming that you have finished the step2 calculations (computed each of the KS-DFT calculations from cp2k)

mkdir all_pdosfiles
for file in $(find ../step2/wd -name '*.pdos'); do cp $file all_pdosfiles/.; echo $file; done

the user only needs to worry about block 1 in the file, excitation_analysis.py
