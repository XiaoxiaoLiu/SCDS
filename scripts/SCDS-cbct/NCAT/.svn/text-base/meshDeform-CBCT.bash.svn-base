#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-6
#$ -j y

NUM=$1;

phaseList=(10 20 30 40 50 60 )
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
type=xyz
mkdir /stage/sharonxx/proj/mskcc/NCAT$NUM/cbct/segmentation
mkdir /stage/sharonxx/proj/mskcc/NCAT$NUM/cbct/segmentation/$type
/home/sharonxx/bin/linux/pDeform  -v  -numPCs 2 -weightShapePrior 0.1  -numIters 25 -stepSizeD 0.3 -gaussSigma 7 -image /stage/sharonxx/proj/mskcc/NCAT$NUM/cbct/image/gray/phase$i.mhd  -pcaStats /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/stats/meshPtsStats-$type/MeshStats.txt  -refMesh /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/stats/meshPtsStats-$type/meanMesh.vtk -output  /stage/sharonxx/proj/mskcc/NCAT$NUM/cbct/segmentation/$type/fit.$i.vtk > /stage/sharonxx/proj/mskcc//NCAT$NUM/cbct/segmentation/$type/fit.$i.vtk.log

echo "done"
