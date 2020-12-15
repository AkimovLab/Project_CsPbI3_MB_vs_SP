## Molecular dynamics calculations

cp2k_md_input_file.inp - this file contains the input for the cp2k md calculations. Actually, this file contains the restart file after the thermalization of the MD. Use this file for directly producing a production MD trajectory.

This files BASIS_MOLOPT, POTENTIAL, and dftd3 contain the Gaussian basis set, pseudopotentials, and detd3 dispersion corrections, respectfully.

The submit file contains bash commands needed to run the cp2k input file on the UB CCR
