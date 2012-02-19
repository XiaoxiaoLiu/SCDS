# bash scripts for running SCDS studies

###################   parmaters setting     #################################
patNo=10146002;  ## patient number: 10146002


echo "Setting parameters..."

## default data storgage settings
binaryFolder=/home/xiaoxiao/work/bin/SCDS/bin
scriptsFolder=/home/xiaoxiao/work/src/SCDS/scripts/SCDS-cbct/diaphragmSCDS

dataFolder=/media/Xiaoxiao_Backup/proj/SCDS/MSKCC_DATA/445_pt$patNo

cbctFolder=$dataFolder/cbct
rcctFolder=$dataFolder/rcct

# output dirs
cbctImageFolder=$dataFolder/cbct/image/original
rcctImageFolder=$dataFolder/rcct/image/original
rcctResampledImageFolder=$dataFolder/rcct/image/original
# rcctAtlasFolder=$rcctFolder/atlas
shapeFolder=$rcctFolder/shape
statsFolder=$rcctFolder/stats
resultsFolder=$dataFolder/results
scdsPredictionFolder=$dataFolder/results/SCDS_predict
scdsWarpFolder=$dataFolder/results/SCDS_warp

mkdir -p $cbctImageFolder
mkdir -p $rcctImageFolder
mkdir -p $shapeFolder
mkdir -p $statsFolder
mkdir -p $resultsFolder
mkdir -p $scdsPredictionFolder
mkdir -p $scdsWarpFolder
mkdir -p $rcctResampledImageFolder




# Converting CBCT Dicom image series *.dcm to Metadata format *.mhd
echo "Converting CBCT Dicom image series *.dcm to Metadata format *.mhd ..."
for i in 00 16 32 50 66 82 ;
do $binaryFolder/DICOMSeriesTo3DImage $cbctFolder/six_phases/cbct1_phase_$i  $cbctImageFolder/phase$i\.mhd ; done


# Converting RCCT Dicom image series *.dcm to Metadata format *.mhd
echo "Converting RCCT Dicom image series *.dcm to Metadata format *.mhd ..."
for i in 0 1 2 3 4 5 6 7 8 9 ; 
do $binaryFolder/DICOMSeriesTo3DImage $rcctFolder/RCCT2_AR/RCCT2_30$i\ANPR $rcctImageFolder/cinePhase$i\0.mhd ; done


#Resampling  RCCT to have the same spacing as CBCT
spacing=(0.908394 0.908394 1.99496)  # CBCT's image sapcing for resampling purpose
for i in 0 1 2 3 4 5 6 7 8 9 ; do $binaryFolder/ResamplingImageBySpacing  3 $rcctImageFolder/cinePhase$i\0.mhd   $rcctResampledImageFolder/cinePhase$i\0.mhd ${spacing[0]} ${spacing[1]} ${spacing[2]}; done

#2. Find out intersection using
$scriptsFolder/intersectionRegionExtract_RCCT_CBCT.m





##############################   Training atlas  ################################################
$scriptsFolder/atlasWerks.bash  Patient$patNo




############################   Predict CBCT deformations (DVFs) ##################################
$scriptsFolder/run_predict_inter.m


#### get motion-corrected average image

# apply the predicted deformation
mkdir  $resultsFolder/SCDS_predict
mkdir  $resultsFolder/SCDS_predict/deformedImage
cd $resultsFolder/SCDS_predict/deformedImage

for n in 00 16 32 50 66 82
do
	txApply -b -i $cbctFolder/image/gray-inter/phase$n.mhd -h $resultsFolder/SCDS_predict/Hfield/pred-hfield-phase$n.mhd  -o deformedPhase$n
done


# for comparison purpose: get the Euclidean mean
AtlasWerks --outputImageFilenamePrefix=averageImage_  --scaleLevel=1 --numberOfIterations=0   *.mhd



### advanced ###
############################  prediction-based atlas formation ##################################

qsub $scriptsFolder/scdsWarp.bash Patient$patNo





###############################        evaluation     ###########################################


evaluateROI.bash

evaluteLungContour.bash ### not done yet





