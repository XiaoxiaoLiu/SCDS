# bash scripts for running jobs on bass.cs.unc.edu
# Xiaoxiao Liu
# Modified  October, 2010



##### Notes about qsub #####
# qsub is a commond used for submitting scripted jobs to bass(.cs.unc.edu), it is used mostly for running multiple jobs
# that are excuted at the same time in different computer nodes on bass.

# For windowns users, the bash scripts in this document could be run directly from cygwin ( or run it on windows "cmd" after minor modifications) except the qsub command.
# Users should take the code out of the input bash script file of each qsub, and strip the actually bash command lines accordingly.  Just take out the first several lines "
# #$ -S /bin/bash #$ -o /home/sharonxx/tmp/   #$ -t 2  #$ -j y  #$ -pe smp"  and replace "$SGE_TASK_ID" with the acutal case numbers for each programs called.





###################   parmaters setting     #################################

patNo=113;  ## patient number: 113

pat$patNo.params.bash  # set threshold, seeds, crop parameters


spacing=(1.52 1.52 1.52)  # CBCT's image sapcing for resampling purpose

## default data storgage settings
scriptsFolder=~sharonxx/scripts/SCDS-cbct/patient
dataFolder=/stage/sharonxx/proj/mskcc/Patient$patNo
cbctFolder=$dataFolder/cbct
cbctImageFolder=$dataFolder/cbct/image

rcctFolder=$dataFolder/rcct
rcctImageFolder=$dataFolder/rcct/image
rcctAtlasFolder=$dataFolder/rcct/atlas
shapeFolder=$dataFolder/rcct/shape
statsFolder=$dataFolder/rcct/stats


resultsFolder=$dataFolder/results
scdsPredictionFolder=$dataFolder/results/SCDS_predict
scdsWarpFolder=$dataFolder/results/SCDS_warp

origDICOMImageFolder=$dataFolder/pt$patNo/cbct1


## Align RRCT with CBCT
#1. Resampling 
for i in 0 1 2 3 4 5 6 7 8 9 ; do ResampleImageBySpacing  3 $rcctImageFolder/gray/cinePhase$i\0.mhd   $rcctImageFolder/gray/cinePhase$i\0.mhd ${spacing[0]} ${spacing[1]} ${spacing[2]}; done

#2. Find out intersection using
$scriptsFolder/intersectionRegionExtract_RCCT_CBCT.m





###################      Binary Segmentation   ################################
#### function binarySeg 
cd $rcctImageFolder
mkdir binary-inter
for i in 00 10 20 30 40 50 60 70 80 90 
do 
	ThresholdImage 3 ./gray-inter/cinePhase$i.mhd ./binary-inter/lung-bin-cinePhase$i.mhd  ${threshold[0]}   ${threshold[1]}
	RollingBall  ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd  0 0  1 1 1 1   25
	ConnectedThreshold ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd   1  2  ${seed[0]}  ${seed[1]}  ${seed[2]} ${seed[3]} ${seed[4]} ${seed[5]}
done

## sanity check the segmentaiton results by looking at the middle slices of all segmentations
#StackSlices StackBinarySegMidSlice.mhd  -1 -1 35  binary-inter/*.mhd



##############  Shape Modeling with particle systems ###################

## 1. Crop the binary image and pad the boundaires
for i in  00 10 20 30 40 50 60 70 80 90 ; 
do CropAndPad  $rcctImageFolder/binary-inter/lung-bin-cinePhase$i.mhd  $shapeFolder/lung-bin-crop-pad-cinePhase$i.mhd  ${crop[0]}  ${crop[1]}  ${crop[2]} ${crop[3]} ${crop[4]} ${crop[5]}  8 8 8 ;
done

# remove small objects
for i in  00 10 20 30 40 50 60 70 80 90; do  BinaryErode  $shapeFolder/lung-bin-crop-pad-cinePhase$i.mhd   $shapeFolder/lung-bin-crop-pad-cinePhase$i.mhd 2; BinaryDilate  $shapeFolder/lung-bin-crop-pad-cinePhase$i.mhd   $shapeFolder/lung-bin-crop-pad-cinePhase$i.mhd 2;done

## 2. Generate distance map for correspondence

qsub $scriptsFolder/distancemap_qsub.bash  $patNo


## 3. Interactively run correspondence 

# write seeds file
$scriptsFolder/writeSeeds.m


mkdir $shapeFolder/xyz
cd  $shapeFolder/xyz

#running with GUI
cp $scriptsFolder/particle.params ./
ShapeWorksShop particle.params 



#running without GUI
cp $scriptsFolder/shapeWorksRun.params ./
ShapeWorksRun ./shapeWorksRun.params




##############  Deformable segmentation on CBCT to extract the surrogate shape ##################

## generate triangle meshes from the original  binary segmenations
qsub $scriptsFolder/vol2surf.bash $patNo 10 # generate the meshes from the original binary volume

## post processing to get the meshes from the correpsonding point sets

cp $scriptsFolder/corresMesh.params  $shapeFolder/xyz
cd $shapeFolder/xyz
ParticleCorrespondencePostprocessingTPS --parameterFileName corresMesh.params 

$scriptsFolder/trainMeshStats.m # train the mesh shape models
mkdir $cbctFolder/segmentation
mkdir $cbctFolder/segmentation/xyz
qsub $scriptsFolder/meshDeform-CBCT.bash $patNo  # segment CBCT





##############################   Training atlas  ################################################
qsub $scriptsFolder/atlasWerks.bash  Patient$patNo




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





