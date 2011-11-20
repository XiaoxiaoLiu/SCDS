#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y


phaseList=(00 10 20 30 40 50 60 70 80 90)
NUM=2;
#running
i=${phaseList[$SGE_TASK_ID-1]}
cd  /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/shape
/home/sharonxx/bin/linux/DistanceMap   lung-bin-crop-pad-cinePhase$i.mhd  lung-bin-cinePhase$i.dis.mhd 15
echo "this is job $ $SGE_TASK_ID"


