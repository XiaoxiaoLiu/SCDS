#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-5
#$ -j y
# registration using fWarp
phaseList=( 10 20 30 50 60 )
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "


echo "this is job  ncat phase$i"
mkdir /stage/sharonxx/proj/mskcc/NCAT/cbct/largeWarp
/home/sharonxx/bin/linux/fWarp  --outputImageFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT/cbct/largeWarp/deformed2p40_$i       --outputHFieldFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT/cbct/largeWarp/hfield-deformed2p40_$i      --scaleLevel=2 --numberOfIterations=250   --scaleLevel=1 --numberOfIterations=100    /stage/sharonxx/proj/mskcc/NCAT/cbct/image/gray/phase40.mhd        /stage/sharonxx/proj/mskcc/NCAT/cbct/image/gray/phase$i.mhd    
echo "done"
