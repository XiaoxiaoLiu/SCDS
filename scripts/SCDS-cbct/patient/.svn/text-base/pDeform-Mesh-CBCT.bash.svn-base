#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-6
#$ -j y


patNo=102
phaseList=(00 16 32 50 66 82)
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "


type=xyz;
/home/sharonxx/bin/linux/pDeform  -v  -numPCs 3 -weightShapePrior 0.1  -numIters 50 -stepSizeD 0.2 -gaussSigma 10 -refMesh  /stage/sharonxx/proj/mskcc/Patient$patNo/rcct/stats/meshPtsStats-$type/meanMesh.vtk  -image /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/image/gray-inter/phase$i.mhd  -pcaStats /stage/sharonxx/proj/mskcc/Patient$patNo/rcct/stats/meshPtsStats-$type/MeshStats.txt -output  /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/segmentation/fit.$i.vtk > /stage/sharonxx/proj/mskcc/Patient$patNo/cbct/segmentation/fit.$i.mesh.log

echo "done"
