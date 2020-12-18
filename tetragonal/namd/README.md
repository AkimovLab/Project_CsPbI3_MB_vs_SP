## namd

In each case of excess electronic excitation energy (here, 0.4 eV), there are 5 folders. The fssh, ida, msdm, and bllz folders contain the scripts for running NA-MD for each of these dynamical cases (ida and msdm are decoherence corrected surface hopping methods and BLLZ is an energy-only (NAC-free) surface hopping method). Each folder contains the script needed to run the dynamics for the particular surface hopping method. The only difference between the `namd.py` script in each folder is a single block of code, Please see the code for more details.

The fit folder contains a script for fitting and plotting the dynamics computed in the fssh, ida, msdm, and bllz folders. It will only work as a black-box if you make no changes in the `namd.py` scripts within the fssh, ida, msdm, and bllz folders. 

As always, 

Be sure to test the first section in each `namd.py` script to make sure that the data is being read properly (are you actually reading in the matricies from step3?).

For video tutorials on using this software, please see [this link](https://github.com/compchem-cybertraining/Tutorials_Libra/blob/master/VIDEOS.md).
