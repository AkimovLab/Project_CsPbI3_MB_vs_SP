## Excitation analysis

This file contains the scripts needed to quantify the degree of configurational mixing for the optimized and thermally averaged systems

In excitation_analysis.py, set "thermal" = False when computing the configurational mixing for the optimized system. Set set "thermal" = True when computing the configurational mixing for the thermally averaged system.

For the thermally averaged system, you will need a folder called "all_logfiles" present in this directy, and it should contain the output (.log) file from each of the LR-TDDFT calculations from cp2k, you can generate this directory by running the following commands (assuming that you have finished the step2 calculations (computed each of the LR-TDDFT calculations from cp2k)

mkdir all_logfiles
for file in $(find ../step2/wd -name '*.log'); do cp $file all_logfiles/.; echo $file; done

the user only needs to worry about block 1 in the file, excitation_analysis.py
