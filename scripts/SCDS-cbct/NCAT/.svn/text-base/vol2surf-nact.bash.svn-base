#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y



phaseList=( 00 10 20 30 40 50 60 70 80 90 )
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "

NUM=2;
echo "this is job  ncat phase$i"
source ~sharonxx/work
# to include the libraries
cd /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/shape
/home/sharonxx/bin/linux/Vol2surf  -input lung-bin-crop-pad-cinePhase$i.mhd  -vtk  -gaussian 15 -label 1  -decimate -targetReduction 85
#98 ---5000 points
echo "done"
