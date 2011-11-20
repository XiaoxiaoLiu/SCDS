#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y



phaseList=( 00 10 20 30 40 50 60 70 80 90 )
phaseNo=$(( $SGE_TASK_ID-1 ))
i=${phaseList[$phaseNo]}
echo "$SGE_TASK_ID;  pahseNO = $phaseNo;  i = $i; "


echo "this is job  ncat phase$i"
source ~sharonxx/work
# to include the libraries
cd /stage/sharonxx/proj/mskcc/NCAT/rcct/shape/old
/home/sharonxx/bin/linux/Vol2surf  -input lung-bin-cinePhase$i.left.mhd  -vtk  -gaussian 32 -label 1 -decimate -targetReduction 2500
echo "done"
