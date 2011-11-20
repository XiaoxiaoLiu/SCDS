patNo=113;

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



# get the frechet mean ROI contour
txApply -b -i $rcctFolder/ROI/cinePhase50.mhd -h $rcctFolder/atlas/hfield_0005.mhd  -o $rcctFolder/ROI/atlas_mean
ThresholdImage 3   $rcctFolder/ROI/atlas_mean.mhd  $rcctFolder/ROI/atlas_mean.mhd 0.5 1




##########################################  reference:  frechet mean of RCCT ####################
mkdir  $scdsWarpFolder/ROI

#SCDS warp OnlyIntensity
mkdir  $scdsWarpFolder/ROI/OnlyIntensity
for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i   $rcctFolder/ROI/atlas_mean.mhd  -h  $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/ROI/OnlyIntensity/phase$i  
ThresholdImage 3 $scdsWarpFolder/ROI/OnlyIntensity/phase$i.mhd$scdsWarpFolder/ROI/OnlyIntensity/phase$i.mhd 0.5 1
done


#SCDS warp with deformation constrains, no shape index, init with Id
mkdir  $scdsWarpFolder/ROI/InitById
for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i   $rcctFolder/ROI/atlas_mean.mhd  -h  $scdsWarpFolder/InitById/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/ROI/InitById/phase$i  
ThresholdImage 3 $scdsWarpFolder/ROI/InitById/phase$i.mhd$scdsWarpFolder/ROI/InitById/phase$i.mhd 0.5 1
done



# SCDS prediction
mkdir $scdsPredictionFolder
mkdir $scdsPredictionFolder/ROI
for i in 00 16 32 50 66 82 
	do txApply -f -i    $cbctFolder/ROI/phase50.mhd   -h $scdsPredictionFolder/Hfield/pred-hfield-phase$i.mhd  -o $scdsPredictionFolder/ROI/phase$i  
ThresholdImage 3  $scdsPredictionFolder/ROI/phase$i.mhd $scdsPredictionFolder/ROI/phase$i.mhd 0.5 1
done



















######################################## reference: phase50 of CBCT #######################
#compute the deformation from phase50 of CBCT






mkdir  $scdsWarpFolder/ROI

#####SCDS warp OnlyIntensity
mkdir  $scdsWarpFolder/ROI/OnlyIntensity

# get the frechet mean ROI contour
txApply -b -i  $cbctFolder/ROI/phase50.mhd -h $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_0003.mhd -o  $cbctFolder/ROI/atlas_mean_Int

for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i    $cbctFolder/ROI/atlas_mean_Int.mhd  -h  $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/ROI/OnlyIntensity/ref-cbct50-phase$i  
ThresholdImage 3 $scdsWarpFolder/ROI/OnlyIntensity/ref-cbct50-phase$i.mhd $scdsWarpFolder/ROI/OnlyIntensity/ref-cbct50-phase$i.mhd 0.5 1
done




######SCDS warp with deformation constrains, no shape index, init with Id
mkdir  $scdsWarpFolder/ROI/InitById

# get the frechet mean ROI contour
txApply -b -i  $cbctFolder/ROI/phase50.mhd -h $scdsWarpFolder/InitById/scdsWarp_Hfield_0003.mhd -o  $cbctFolder/ROI/atlas_mean_InitById

for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i  $cbctFolder/ROI/atlas_mean_InitById.mhd  -h  $scdsWarpFolder/InitById/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/ROI/InitById/ref-cbct50-phase$i  
ThresholdImage 3 $scdsWarpFolder/ROI/InitById/ref-cbct50-phase$i.mhd  $scdsWarpFolder/ROI/InitById/ref-cbct50-phase$i.mhd 0.5 1
done



##### SCDS prediction
mkdir $scdsPredictionFolder
mkdir $scdsPredictionFolder/ROI

# get the frechet mean ROI contour
txApply -b -i  $cbctFolder/ROI/phase50.mhd -h $scdsPredictionFolder/Hfield/pred-hfield-phase50.mhd -o  $cbctFolder/ROI/atlas_mean_pred

for i in 00 16 32 50 66 82 
do txApply -f -i    $cbctFolder/ROI/atlas_mean_pred.mhd  -h $scdsPredictionFolder/Hfield/pred-hfield-phase$i.mhd  -o $scdsPredictionFolder/ROI/ref-cbct50-phase$i  
ThresholdImage 3  $scdsPredictionFolder/ROI/ref-cbct50-phase$i.mhd $scdsPredictionFolder/ROI/ref-cbct50-phase$i.mhd 0.5 1
done


