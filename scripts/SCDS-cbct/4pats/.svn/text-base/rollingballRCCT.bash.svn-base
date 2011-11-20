#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y

patList=( 21 )
# 20 25 31 
phaseList=( 00 10 20 30 40 50 60 70 80 90 )
patNo=$(( ($SGE_TASK_ID-1)/10))
phaseNo=$(( $SGE_TASK_ID-$patNo*10-1 ))
i=${phaseList[$phaseNo]}
j=${patList[$patNo]}

echo "this is job pat$j phase$i"

cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/image

/home/sharonxx/bin/linux/RollingBall  ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-RB-cinePhase$i.mhd  0 0   1 1 1  0   25

