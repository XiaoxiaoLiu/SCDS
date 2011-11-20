#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-10
#$ -j y


patList=( 21)
#20 21 25 31 ) 
phaseList=( 00 10 20 30 40 50 60 70 80 90 )
patNo=$(( ($SGE_TASK_ID-1)/10))
phaseNo=$(( $SGE_TASK_ID-$patNo*10-1 ))
i=${phaseList[$phaseNo]}
j=${patList[$patNo]}
echo "$SGE_TASK_ID; patNO = $patNo ; pahseNO = $phaseNo;  i = $i; j= $j"
echo "this is job pat$j phase$i"
echo "this is job pat$j phase$i"
mkdir /stage/sharonxx/proj/mskcc/Patient$j/rcct/shape
/home/sharonxx/bin/linux/DistanceMap /stage/sharonxx/proj/mskcc/Patient$j/rcct/shape/lung-bin-cinePhase$i.mhd  /stage/sharonxx/proj/mskcc/Patient$j/rcct/shape/lung-bin-cinePhase$i.DT.mhd 9

echo "this is job $ $SGE_TASK_ID"


