while read trait
do
    mkdir -p farmCPU_${trait}
    echo ${trait} > farmCPU_${trait}/params.txt

    cd farmCPU_${trait}
    nohup Rscript ../run_farmCPU.R 5 > farmCPU_${trait}_5PCs.log &
    cd ../
done < traits.txt



