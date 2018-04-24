# GWAS_WiDiv-942_farmCPU

Contains scripts for GWAS using the farmCPU software described in Liu et al. (2016; doi:10.1371/journal.pgen.1005767).  The WiDiv-942 is a diversity panel described in Hirsch et al. (2014; doi:10.1105/tpc.113.119982) and more recently in Mazaheri et al. (2018; _paper doi_).  The imputed genotypic data for the WiDiv-942 is available at doi:10.5061/dryad.n0m260p.

* **1_prep_dirs_run_farmCPU.sh**: this script creates a results directory for each trait, then calls run_farmCPU.R and logs the output.
* **run_farmCPU.R**: this script runs farmCPU.  It accepts a single number as an argument, specifying the number of principal components to use as covariates.  Can also run permutations to determine the entry threshold for adding SNPs to the model during the first iteration.
