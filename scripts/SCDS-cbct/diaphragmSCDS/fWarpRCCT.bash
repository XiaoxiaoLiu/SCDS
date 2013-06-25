#!/bin/bash
# registration using fWarp

################### Setting parmaters  #################################

echo "Setting parameters..."

patNo=Pt27;  ## patient number: 10146002

## default data storgage settings

binaryFolder=/home/xiaoxiao/work/bin/SCDS/bin
scriptsFolder=/home/xiaoxiao/work/src/SCDS/scripts/SCDS-cbct/diaphragmSCDS

dataFolder=/media/Xiaoxiao_Backup/proj/SCDS/MSKCC_DATA/$patNo

cbctFolder=$dataFolder/cbct
rcctFolder=$dataFolder/rcct

rcctImageFolder=$dataFolder/rcct/image
echo "Setting parameters done!"
echo
echo " Training atlas RCCT ..."


##############    Training atlas RCCT      ###################
mkdir $rcctFolder/largeWarp
for i in 05 15 25 35 45 65 75 85 95
do $binaryFolder/fWarp  \
--outputImageFilenamePrefix=$rcctFolder/largeWarp/deformed_cinePhase$i \
--outputHFieldFilenamePrefix=$rcctFolder/largeWarp/hfield-cinePhase$i \
--intensityWindowMin=-1024   --intensityWindowMax=1024   \
--scaleLevel=4 --numberOfIterations=300 \
--alpha=0.01  --beta=0.01  --gamma=0.001  \
--scaleLevel=2 --numberOfIterations=200 \
--alpha=0.01  --beta=0.01  --gamma=0.001  \
--scaleLevel=1 --numberOfIterations=100 \
--alpha=0.01  --beta=0.01   --gamma=0.001  \
$rcctImageFolder/inter/cinePhase55.mhd  \
$rcctImageFolder/inter/cinePhase$i.mhd  \
> $rcctFolder/largeWarp/fWarp$i.log 2>&1
done

echo "done"
