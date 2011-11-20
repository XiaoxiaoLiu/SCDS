#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-6
#$ -j y



phaseList=(10 20 30 40 50 60 )
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "


echo "this is job  ncat phase$i"
#type=2curvature;
type=xyz;
/home/sharonxx/bin/linux/pDeform  -v  -numPCs 3 -weightShapePrior 10  -numIters 15 -stepSizeD 0.4 -gaussSigma 10 -image /stage/sharonxx/proj/mskcc/NCAT/cbct/image/gray/phase$i.mhd  -pcaStats /stage/sharonxx/proj/mskcc/NCAT/rcct/stats/PtsStats-$type/PtsStats_NCAT.txt -output  /stage/sharonxx/proj/mskcc/NCAT/cbct/segmentation/$type/fit.$i.lpts > /stage/sharonxx/proj/mskcc/NCAT/cbct/segmentation/$type/fit.$i.pts.log

echo "done"
