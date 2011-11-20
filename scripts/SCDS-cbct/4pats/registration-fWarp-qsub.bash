#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-27
#$ -j y
# registration using fWarp
patList=( 20 21 25 31 ) 
phaseList=( 00 10 20 30 40 60 70 80 90 )
patNo=$(( ($SGE_TASK_ID-1)/9))
phaseNo=$(( $SGE_TASK_ID-$patNo*9-1 ))
i=${phaseList[$phaseNo]}
j=${patList[$patNo]}
echo "$SGE_TASK_ID; patNO = $patNo ; pahseNO = $phaseNo;  i = $i; j= $j"


echo "this is job pat$j phase$i"
mkdir /stage/sharonxx/proj/mskcc/Patient$j/largeWarp
/home/sharonxx/bin/linux/fWarp  --outputImageFilenamePrefix=/stage/sharonxx/proj/mskcc/Patient$j/largeWarp/deformed_cinePhase$i       --outputHFieldFilenamePrefix=/stage/sharonxx/proj/mskcc/Patient$j/largeWarp/hfield-cinePhase$i                                      --scaleLevel=4 --numberOfIterations=200                                         --scaleLevel=2 --numberOfIterations=50                                          --scaleLevel=1 --numberOfIterations=20                                            /stage/sharonxx/proj/mskcc/Patient$j/rcct/image/gray-inter/cinePhase50.mhd  /stage/sharonxx/proj/mskcc/Patient$j/rcct/image/gray-inter/cinePhase$i.mhd
echo "done"
