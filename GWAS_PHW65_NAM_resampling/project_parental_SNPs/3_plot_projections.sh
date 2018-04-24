for POP in PHN11 MoG Mo44; do
    for CHR in $( seq 10 ); do
	nohup Rscript plot_projections_and_reformat.R $POP $CHR &> ${POP}_chr${CHR}_plotAndReformat.log &
    done
done
