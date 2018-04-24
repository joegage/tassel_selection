# selection

Contains scripts used to calculate Fst, XP-EHH, and XP-CLR between inbreds derived from cycle 0 of the BSSS population and Stiff Stalk ex-PVPs.

## XPCLR

* **1_create_files.R**: pulls the appropriate individuals out of the WiDiv-942 genotype file and formats them for XP-CLR calculation.
* **2_run_xpclr.sh**: calculates XP-CLR along each chromosome.  This is done using a modified version of the `XPCLR` program created by Chen et al. (2010, doi:10.1101/gr.100545.109).  The modification was done by Hufford et al. (2012, doi:10.1038/ng.2309) to account for windows with no SNP data.
* **3_eval_SNP_num.R**: creates plots to evaluate the number of SNPs per window on each chromosome.  Also plots XP-CLR values against physical position on each chromosome.
* **4_bin_XPCLR.R**: calculates maximum XP-CLR scores in 10kb windows, writes the results to file.

## XPEHH_Fst

* **1_mod_files.R**: modifies the XP-CLR genotype and map files to be suited for running in `hapbin`.
* **2_run_xpehh.sh**: calculates XP-EHH for each chromosome using `hapbin` (Mclean et al. 2015, doi:10.1093/molbev/msv172).
* **3_fst.R**: calculates Fst at each SNP along each chromosome using the `vectorFst` function (Beissinger et al. 2014, doi:10.1534/genetics.113.160655, available at http://beissingerlab.github.io/docs/vectorFst.R)
* **4_bin_XPEHH_Fst.R**: calculates maximum XP-EHH and Fst in 10kb windows, writes the results to file.
