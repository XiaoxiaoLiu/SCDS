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

cbctImageFolder=$cbctFolder/image
echo "Setting parameters done!"
echo
echo " Training atlas RCCT ..."


##############    Training atlas RCCT      ###################
mkdir $cbctFolder/largeWarp

for i in 00 16 32 66 82
  do $binaryFolder/fWarp  \
--outputImageFilenamePrefix=$cbctFolder/largeWarp/deformed_cinePhase$i \
--outputHFieldFilenamePrefix=$cbctFolder/largeWarp/hfield-cinePhase$i \
--scaleLevel=4 --numberOfIterations=200 \
--scaleLevel=2 --numberOfIterations=100 \
--scaleLevel=1 --numberOfIterations=50 \
$cbctImageFolder/inter/phase50.mhd  \
$cbctImageFolder/inter/phase$i.mhd  \
> $cbctFolder/largeWarp/fWarp$i.log 2>&1

  done

echo "done"
