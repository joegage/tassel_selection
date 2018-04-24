# GWAS_PHW65_NAM_resampling
Contains scripts used for resampling GWAS analysis of the PHW65 NAM population.  This follows the mapping methods detailed in Tian et al. (2011; doi:10.1038/ng.746).  The PHW65 NAM is a nested mapping population consisting of three biparental populations of 200 individuals each.  The reference line (common parent) is PHW65, while the three founder lines are PHN11, Mo44, and MoG.

* **1_run_StepwiseAdditiveModelFitterPlugin.sh**: runs the first step of NAM analysis by performing stepwise regression on approximately 9k GBS markers.  Chromosomeare-specific residuals are written to individual files for use in the resampling GWAS step.
* **2_run_ResamplingGWASPlugin.sh**: runs the second step of NAM analysis.  Uses the chromosome-specific residuals from the first step as response variables, subsamples 80% of the individuals in each population, and adds SNPs to a forward regression model. This process is repeated 100 times and the SNPs added at each step of each subsampling iteration are written to file.
* **3_permute_resamplingGWAS.sh**: runs permutations in order to determine the RMIP threshold that can be considered significant.  The same model as in step 2 (above) is run, but with the chromosome-specific residuals randomized.  An appropriate RMIP threshold for the true analysis can then be calculated based on the RMIPs of SNPs from the permutations.

**project_parental_SNPs**: this directory contains scripts used to project high-density parental SNPs onto the progeny, using the ~9k GBS markers as a scaffold.
