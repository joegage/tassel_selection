## XPCLR=~/applications/XPCLR/bin/XPCLR # BAD software with bug
XPCLR=~/applications/src_NA/XPCLR

# ORIGINAL:
# gWIN=0.001
# snpWIN=50
# GRID=1000
# COR=0.7

gWIN=0.001 ## 0.001
snpWIN=50
GRID=100 ## 100
COR=0.7


for CHROM in 1 2 3 4 5 6 7 8 9 10; do
	
	BSSS=BSSS_chrom${CHROM}.txt
	PVP=exPVP_chrom${CHROM}.txt
	SNP=snpInfo_chrom${CHROM}.txt
	
	OUTFILE=xpclr_chrom${CHROM}
	
	nohup $XPCLR -xpclr $PVP $BSSS $SNP $OUTFILE -w1 $gWIN $snpWIN $GRID $CHROM -p1 $COR &> xpclr_chrom${CHROM}.log &

	if [ $CHROM = 5 ]; then
	   wait
	fi
	
done
