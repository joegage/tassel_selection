for POP in PHN11 MoG Mo44; do
    for CHR in $( seq 10 ); do
	nohup Rscript project_parental_SNPs.R $POP $CHR &> ${POP}_chr${CHR}_projected10M.log &
    done
done
