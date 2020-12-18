## step3

This file contains the scripts needed convert overlap, time-overlap, and energies in step2 (which are in the KS basis) to the single Slater determinant (single particle (SP)) basis as the many determinant (many-body (MB)) bases.

For the file `build_bases.py` to run properly, you will need a folder called `all_logfiles` present in the `excitation_analysis` directory, and it should contain the output (`.log`) file from each of the LR-TDDFT calculations from CP2K. Please see the `README.md` in that folder for further details. 

Below are a couple tips for running `build_bases.py`.

Make sure the the variable `data_dim` is set to the total number of rows in the matricies formed in step2. That is, the number of rows in the files that contain the overlap, time-overlap, and energies in step2. Also make sure that the values such as: `min_band_ks`, `max_band_ks`, and `ks_homo_index` are defined correctly. These values are the same ones you defined in `step2/submit_template.slm`.

Always test the first section to make sure that the data is being read properly (are you actually reading in the matricies from step2?).

For video tutorials on using this software, please see [this link](https://github.com/compchem-cybertraining/Tutorials_Libra/blob/master/VIDEOS.md).
