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
/home/sharonxx/bin/linux/pDeform -v  -numPCs 5 -weightShapePrior 10  -numIters 15 -stepSizeD 0.4 -gaussSigma 10 -image /stage/sharonxx/proj/mskcc/NCAT/cbct/image/gray/phase$i.mhd  -pcaStats /stage/sharonxx/proj/mskcc/NCAT/rcct/stats/meshPtsStats-$type/meshPtsStats_NCAT.txt -refMesh /stage/sharonxx/proj/mskcc/NCAT/rcct/stats/meshPtsStats-$type/meanMesh_NCAT.vtk -output  /stage/sharonxx/proj/mskcc/NCAT/cbct/segmentation/$type/seg.$i.vtk > /stage/sharonxx/proj/mskcc/NCAT/cbct/segmentation/$type/seg.$i.vtk.log

echo "done"
