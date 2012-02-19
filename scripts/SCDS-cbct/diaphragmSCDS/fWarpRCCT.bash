#!/bin/bash
# registration using fWarp

################### Setting parmaters  #################################

echo "Setting parameters..."

patNo=10146002;  ## patient number: 10146002

## default data storgage settings

binaryFolder=/home/xiaoxiao/work/bin/SCDS/bin
scriptsFolder=/home/xiaoxiao/work/src/SCDS/scripts/SCDS-cbct/diaphragmSCDS

dataFolder=/media/Xiaoxiao_Backup/proj/SCDS/MSKCC_DATA/445_pt$patNo

cbctFolder=$dataFolder/cbct
rcctFolder=$dataFolder/rcct

rcctImageFolder=$dataFolder/rcct/image
echo "Setting parameters done!"
echo
echo " Training atlas RCCT ..."


##############    Training atlas RCCT      ###################
mkdir $rcctFolder/largeWarp

for i in 00 10 20 30 40 60 70 80 90
  do $binaryFolder/fWarp  \
--outputImageFilenamePrefix=$rcctFolder/largeWarp/deformed_cinePhase$i \
--outputHFieldFilenamePrefix=$rcctFolder/largeWarp/hfield-cinePhase$i \
--scaleLevel=4 --numberOfIterations=200 \
--scaleLevel=2 --numberOfIterations=100 \
--scaleLevel=1 --numberOfIterations=100 \
$rcctImageFolder/inter/cinePhase50.mhd  \
$rcctImageFolder/inter/cinePhase$i.mhd  \
> $rcctFolder/largeWarp/fWarp$i.log 2>&1

  done

echo "done"
