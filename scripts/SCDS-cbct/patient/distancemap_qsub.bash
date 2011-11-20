#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y

patNo=$1;
phaseList=(00 10 20 30 40 50 60 70 80 90)
#running
i=${phaseList[$SGE_TASK_ID-1]}
cd  /stage/sharonxx/proj/mskcc/Patient$patNo/rcct/shape
/home/sharonxx/bin/linux/DistanceMap   lung-bin-crop-pad-cinePhase$i.mhd  lung-bin-cinePhase$i.dis.mhd 15
echo "this is job $ $SGE_TASK_ID"


