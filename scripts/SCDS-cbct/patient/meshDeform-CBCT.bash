#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-6
#$ -j y

patNo=$1
phaseList=(00 16 32 50 66 82)
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "


type=xyz;
mkdir /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/segmentation;
/home/sharonxx/bin/linux/pDeform  -v  -numPCs 2 -weightShapePrior 1.0  -numIters 20 -stepSizeD 0.4 -gaussSigma 5 -image /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/image/gray-inter/phase$i.mhd  -pcaStats /stage/sharonxx/proj/mskcc/Patient$patNo/rcct/stats/meshPtsStats-$type/MeshStats.txt  -refMesh /stage/sharonxx/proj/mskcc/Patient$patNo/rcct/stats/meshPtsStats-$type/meanMesh.vtk -output  /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/segmentation/$type/fit.$i.vtk > /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/segmentation/$type/fit.$i.vtk.log

echo "done"
