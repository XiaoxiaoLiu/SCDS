#!/bin/bash

# bash scripts for running SCDS studies at MSKCC
# Author: Oleksandr Dzyubak
# Date: January 05, 2012

# Modified: Oleksandr Dzyubak
# Last Modification:  January 05, 2012

## default data storgage settings

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

mkdir  $rcctFolder/atlas
# cd     $rcctFolder/atlas
rcctAtlasFolder=$rcctFolder/atlas

$binaryFolder/AtlasWerks \
--outputImageFilenamePrefix=$rcctAtlasFolder/averageImage_  \
--outputDeformedImageFilenamePrefix=$rcctAtlasFolder/deformedImage_  \
--outputHFieldFilenamePrefix=$rcctAtlasFolder/hfield_   \
--scaleLevel=4  \
--numberOfIterations=1 \
--scaleLevel=2  \
--numberOfIterations=1 \
--scaleLevel=1  \
--numberOfIterations=1 \
 $rcctImageFolder/inter/cinePhase*0.mhd > $rcctAtlasFolder/atlasWerks.log 2>&1

echo
echo "Processing images complete!"
