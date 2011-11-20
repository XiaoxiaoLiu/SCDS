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
mkdir /stage/sharonxx/proj/mskcc/NCAT3/cbct/segmentation/$type
/home/sharonxx/bin/linux/pDeform  -v  -numPCs 2 -weightShapePrior 1  -numIters 15 -stepSizeD 0.4 -gaussSigma 10 -image /stage/sharonxx/proj/mskcc/NCAT3/cbct/image/gray/phase$i.mhd  -pcaStats /stage/sharonxx/proj/mskcc/NCAT1/rcct/stats/meshPtsStats-$type/MeshStats.txt    -refMesh /stage/sharonxx/proj/mskcc/NCAT2/rcct/stats/meshPtsStats-xyz/meanMesh.vtk -output  /stage/sharonxx/proj/mskcc/NCAT3/cbct/segmentation/$type/fit.$i.NCAT1.vtk > /stage/sharonxx/proj/mskcc/NCAT3/cbct/segmentation/$type/fit.$i.NCAT1.vtk.log

echo "done"

