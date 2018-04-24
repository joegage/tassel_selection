# selection

Contains scripts used to calculate Fst, XP-EHH, and XP-CLR between inbreds derived from cycle 0 of the BSSS population and Stiff Stalk ex-PVPs.

## XPEHH_Fst

* **1_mod_files.R**: modifies the XP-CLR genotype and map files to be suited for running in `hapbin`.
* **2_run_xpehh.sh**: calculates XP-EHH for each chromosome using `hapbin` (Mclean et al. 2015, doi:10.1093/molbev/msv172).
* **3_fst.R**: calculates Fst at each SNP along each chromosome using the `vectorFst` function (Beissinger et al. 2014, doi:10.1534/genetics.113.160655, available at http://beissingerlab.github.io/docs/vectorFst.R)
* **4_bin_XPEHH_Fst.R**: calculates maximum XP-EHH and Fst in 10kb windows, writes the results to file.
