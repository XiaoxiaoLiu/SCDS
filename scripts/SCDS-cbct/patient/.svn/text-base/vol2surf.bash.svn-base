#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y


patNo=$1;
guassVar=$2;

phaseList=( 00 10 20 30 40 50 60 70 80 90 )
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "

echo "this is job  ncat phase$i"

# to include the libraries
export LD_LIBRARY_PATH=/home/sharonxx/lib/i686

cd /stage/sharonxx/proj/mskcc/Patient$patNo/rcct/shape
/home/sharonxx/bin/linux/Vol2surf  -input lung-bin-crop-pad-cinePhase$i.mhd  -vtk  -gaussian $guassVar -label 1
# -decimate -targetReduction 50
#98 ---5000 points
echo "done"
