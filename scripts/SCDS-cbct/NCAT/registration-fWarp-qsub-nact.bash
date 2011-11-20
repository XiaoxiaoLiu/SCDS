#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-6
#$ -j y
# registration using fWarp
phaseList=( 00 10 20 30 40 60 70 80 90 )
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "

NUM=3;
echo "this is job  ncat phase$i"
mkdir /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/largeWarp
/home/sharonxx/bin/linux/fWarp  --outputImageFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/largeWarp/deformed_to_p50_cinePhase$i       --outputHFieldFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/largeWarp/hfield-deformed_to_p50-cinePhase$i   --writeInverseHFields=TRUE   --scaleLevel=4  --alpha=0.1 --beta=0.1  --gamma=0.01  --numberOfIterations=1  --scaleLevel=2 --alpha=0.01 --beta=0.01  --gamma=0.001 --numberOfIterations=50    --scaleLevel=1 --alpha=0.01 --beta=0.01  --gamma=0.001 --numberOfIterations=25     /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/image/gray-inter/cinePhase50.mhd        /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/image/gray-inter/cinePhase$i.mhd  




#/home/sharonxx/bin/linux/fWarp  --outputImageFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT/rcct/largeWarp/inv-deformed_cinePhase$i       --outputHFieldFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT/rcct/largeWarp/inv-hfield-cinePhase$i                                      --scaleLevel=4 --numberOfIterations=200                                         --scaleLevel=2 --numberOfIterations=100                                          --scaleLevel=1 --numberOfIterations=50                                           /stage/sharonxx/proj/mskcc/NCAT/rcct/image/gray-inter/cinePhase$i.mhd  /stage/sharonxx/proj/mskcc/NCAT/rcct/image/gray-inter/cinePhase50.mhd  
echo "done"
