XPEHH=~/applications/hapbin/build/xpehhbin

for CHR in $( seq 10 ); do

    echo $CHR
    
    BSSS=BSSS_chrom${CHR}_xpehh.txt
    PVP=exPVP_chrom${CHR}_xpehh.txt
    INFO=snpInfo_chrom${CHR}_xpehh.txt
    OUT=xpehh_chrom${CHR}

    nohup $XPEHH --hapA $BSSS --hapB $PVP --map $INFO --minmaf 0 --scale 20000 --out $OUT.txt &> $OUT.log &
    
     wait

done
