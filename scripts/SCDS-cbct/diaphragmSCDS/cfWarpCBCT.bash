#!/bin/bash

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

##############    Training atlas RCCT      ###################
#mkdir $dataFolder/results/cfWarp

for i in 00 16 32 66 82
  do $binaryFolder/cfWarp  \
--outputImageFilenamePrefix=$dataFolder/results/cfWarp/deformed_cinePhase$i \
--outputHFieldFilenamePrefix=$dataFolder/results/cfWarp/hfield-cinePhase$i \
--inputHFieldFile=$dataFolder/results/Diaphgram_predict/Hfield/pred-hfield-phase$i.mhd \
--HWeights 1.0  1.0  1.0 \
--scaleLevel=4 --numberOfIterations=100 \
--scaleLevel=2 --numberOfIterations=50 \
--scaleLevel=1 --numberOfIterations=50 \
$cbctFolder/image/inter/phase50.mhd  \
$cbctFolder/image/inter/phase$i.mhd  \
> $dataFolder/results/cfWarp/fWarp$i.log 2>&1
  done

echo "done"
