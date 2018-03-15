
nPerm=5 # Number of permutations to run
limit=1.0E-4 # Entry limit for SNPs to be added to the resampling model

for PERM in $( seq $nPerm ); do

    for CHR in $( seq 10 ); do

	genoFile=../project_parental_SNPs/allPops_chr${CHR}_projected10M_TASSEL.txt
	importString="-fork1 -importGuess $genoFile"
	resampleString=""
	runForkString="-runfork1"
	
	for NUM in $( seq $( wc -l traits.txt | cut -f1 -d" " )); do
	    
	    TRAIT=$( head -$NUM traits.txt | tail -1 )
	    NUMplusOne=$( expr $NUM + 1 )

	    mkdir -p ./permutations/$TRAIT

	    outBase=./permutations/${TRAIT}/permute_${TRAIT}_chr${CHR}_perm${PERM}
	    residFile=./results/${TRAIT}/stepwise_${TRAIT}_chr_${CHR}.txt
	    permPhenoFile=./permutations/${TRAIT}/shuffleResid_chr${CHR}_perm${PERM}.txt

	    # Make permuted phenotype file
	    head -3 $residFile > $permPhenoFile
	    tail -n+4 $residFile | cut -f1-2 > tmp.txt
	    tail -n+4 $residFile | cut -f3 | shuf | paste tmp.txt - >> $permPhenoFile

	    # Construct import section of TASSEL command
	    importString=$( echo $importString -fork${NUMplusOne} -r ${permPhenoFile} )
	    # Construct portion of TASSEL command to merge datasets and run the plugin
	    resampleString=$( echo $resampleString -combine -input1 -input${NUMplusOne} -intersect -ResamplingGWASPlugin -enterLimit $limit -endPlugin -export $outBase.txt)
	    # Construct portion of TASSEL command to run each fork
	    runForkString=$( echo $runForkString -runfork$NUMplusOne )
	    
	done

	command="./TASSEL5/run_custom_pipeline.pl -Xms50g -Xmx100g -debug $importString $resampleString $runForkString"
	nohup $command &> permutations/permute_chr${CHR}_perm${PERM}.log &

	# Only run half of the chromosomes at 
	if [ $CHR = 5 ]
	then
	    wait
	fi

    done

    # Only run one permutation at a time
    wait
    
done
