#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-6
#$ -j y


patNo=105
phaseList=(00 16 32 50 66 82)
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "


type=xyz;
mkdir /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/segmentation;
/home/sharonxx/bin/linux/pDeform  -v  -numPCs 3 -weightShapePrior 0.1  -numIters 25 -stepSizeD 0.2 -gaussSigma 7 -image /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/image/gray-inter/phase$i.mhd  -pcaStats /stage/sharonxx/proj/mskcc/Patient$patNo/rcct/stats/PtsStats-$type/PtsStats.txt -output  /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/segmentation/$type/fit.$i.lpts > /stage/sharonxx/proj/mskcc//Patient$patNo/cbct/segmentation/$type/fit.$i.pts.log

echo "done"
