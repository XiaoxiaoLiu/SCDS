#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y
#$ -pe smp 2
# registration using fWarp

DATANAME=test1
phaseList=(01 02 03 04 05 06 07 08 09 10)
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "

echo "this is job  ncat phase$i"
mkdir /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/largeWarp
/home/sharonxx/bin/linux/fWarp  --outputImageFilenamePrefix=/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/largeWarp/deformed2p5_$i       --outputHFieldFilenamePrefix=/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/largeWarp/hfield-deformed2p5_$i  --scaleLevel=2  --numberOfIterations=200 --scaleLevel=1  --numberOfIterations=100  /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/orig/sphere05.mhd   /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/orig/sphere$i.mhd    > /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/largeWarp/fWarp_deformed2p5_$i.log 2>&1 
echo "done"
