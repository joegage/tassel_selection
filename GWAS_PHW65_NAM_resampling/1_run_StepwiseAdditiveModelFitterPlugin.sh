while read TRAIT; do
   echo $TRAIT
   
    mkdir -p ./results/$TRAIT
    rm -f ./results/$TRAIT/*

    phenoFile=phenotypes/blups_${TRAIT}_forTASSEL.txt
    genoFile=imputed_parents_ChrAll_AllPops_withSampleNames_FINAL.hmp.txt 

    stepwiseCommand="./TASSEL5/run_custom_pipeline.pl -Xmx10g -debug \
    	-fork1 \
    	-h $genoFile \
    	-fork2 \
    	-importGuess $phenoFile \
    	-combine3 -input1 -input2 -intersect \
    	-StepwiseAdditiveModelFitterPlugin \
    	  -criterion pval \
    	  -isNested true \
    	  -usePerm true \
    	  -nPerm 1000 \
    	  -residuals true \
    	  -effectsPrescan true \
    	  -saveToFile true \
    	  -savePath ./results/$TRAIT/stepwise \
    	-endPlugin \
    	-runfork1 -runfork2"

    $stepwiseCommand > results/${TRAIT}/stepwise_${TRAIT}.log

done < traits.txt
