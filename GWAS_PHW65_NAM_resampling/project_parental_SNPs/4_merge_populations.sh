for CHR in $( seq 10 ); do

    cat PHN11_chr${CHR}_projected10M_TASSEL.txt > allPops_chr${CHR}_projected10M_TASSEL.txt
    tail -n +3 Mo44_chr${CHR}_projected10M_TASSEL.txt >> allPops_chr${CHR}_projected10M_TASSEL.txt
    tail -n +3 MoG_chr${CHR}_projected10M_TASSEL.txt >> allPops_chr${CHR}_projected10M_TASSEL.txt

done

    
