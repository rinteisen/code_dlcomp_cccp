mcc -R -nojvm -R -nodisplay -R -nosplash -R -nodesktop -mv ./main_v4_cccp_sqp.m

for ch in {1..10}; do
    for i in {100,200,400}; do
    	#statements
        echo "
        ### Job names
        #PBS -N DL_cccp-$ch-$p
        ### out files
        #PBS -e ./log/DL_cccp-$ch-$p.err
        #PBS -o ./log/DL_cccp-$ch-$p.log
        ### put the job to which queue (qwork)
        #PBS -q qwork" > ./DL_cccp.sh
        echo ' 
        echo Working directory is $PBS_O_WORKDIR
        cd $PBS_O_WORKDIR
        echo noCoMPning on host `hostname`
        echo Start time is `date`
        time1=`date +%s`
        echo Directory is `pwd`' >> ./DL_cccp.sh
        echo " 
        ./main_v4_cccp_sqp $ch 20 6" >> ./DL_cccp.sh
        echo '
        echo End time is `date`
        time2=`date +%s`
        echo Computing time is `echo $time2-$time1 | bc` sec
        ' >> ./DL_cccp.sh
        qsub DL_cccp.sh
    done
done