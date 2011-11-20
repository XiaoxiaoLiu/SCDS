

####################################### not working yet################################
##################### CNT quality is bad ############################

# convert CNT to GTV binary image
# align the GTV to the intersected image
$scriptsFolder/writeCNT2Meta.m

patNo=105;

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



# get the frechet mean CNT contour
txApply -b -i $rcctFolder/CNT/cinePhase50.mhd -h $rcctFolder/atlas/hfield_0005.mhd  -o $rcctFolder/CNT/atlas_mean
ThresholdImage 3   $rcctFolder/CNT/atlas_mean.mhd  $rcctFolder/CNT/atlas_mean.mhd 0.5 1



#######################################  not in use##############################
##########################################  reference:  frechet mean of RCCT ####################
mkdir  $scdsWarpFolder/CNT

#SCDS warp OnlyIntensity
mkdir  $scdsWarpFolder/CNT/OnlyIntensity
for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i   $rcctFolder/CNT/atlas_mean.mhd  -h  $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/CNT/OnlyIntensity/phase$i  
ThresholdImage 3 $scdsWarpFolder/CNT/OnlyIntensity/phase$i.mhd$scdsWarpFolder/CNT/OnlyIntensity/phase$i.mhd 0.5 1
done


#SCDS warp with deformation constrains, no shape index, init with Id
mkdir  $scdsWarpFolder/CNT/InitById
for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i   $rcctFolder/CNT/atlas_mean.mhd  -h  $scdsWarpFolder/InitById/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/CNT/InitById/phase$i  
ThresholdImage 3 $scdsWarpFolder/CNT/InitById/phase$i.mhd$scdsWarpFolder/CNT/InitById/phase$i.mhd 0.5 1
done



# SCDS prediction
mkdir $scdsPredictionFolder
mkdir $scdsPredictionFolder/CNT
for i in 00 16 32 50 66 82 
	do txApply -f -i    $cbctFolder/CNT/phase50.mhd   -h $scdsPredictionFolder/Hfield/pred-hfield-phase$i.mhd  -o $scdsPredictionFolder/CNT/phase$i  
ThresholdImage 3  $scdsPredictionFolder/CNT/phase$i.mhd $scdsPredictionFolder/CNT/phase$i.mhd 0.5 1
done

#######################################  not in use##############################

















######################################## reference: phase50 of CBCT #######################
#compute the deformation from phase50 of CBCT




##### SCDS prediction
mkdir $scdsPredictionFolder
mkdir $scdsPredictionFolder/CNT

# get the frechet mean CNT contour
txApply -b -i  $cbctFolder/CNT/phase50.mhd -h $scdsPredictionFolder/Hfield/pred-hfield-phase50.mhd -o  $cbctFolder/CNT/atlas_mean_pred

for i in 00 16 32 50 66 82 
do txApply -f -i    $cbctFolder/CNT/atlas_mean_pred.mhd  -h $scdsPredictionFolder/Hfield/pred-hfield-phase$i.mhd  -o $scdsPredictionFolder/CNT/ref-cbct50-phase$i  
ThresholdImage 3  $scdsPredictionFolder/CNT/ref-cbct50-phase$i.mhd $scdsPredictionFolder/CNT/ref-cbct50-phase$i.mhd 0.5 1
done





mkdir  $scdsWarpFolder/CNT

#####SCDS warp OnlyIntensity
mkdir  $scdsWarpFolder/CNT/OnlyIntensity

# get the frechet mean CNT contour
txApply -b -i  $cbctFolder/CNT/phase50.mhd -h $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_0003.mhd -o  $cbctFolder/CNT/atlas_mean_Int

for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i    $cbctFolder/CNT/atlas_mean_Int.mhd  -h  $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/CNT/OnlyIntensity/ref-cbct50-phase$i  
ThresholdImage 3 $scdsWarpFolder/CNT/OnlyIntensity/ref-cbct50-phase$i.mhd $scdsWarpFolder/CNT/OnlyIntensity/ref-cbct50-phase$i.mhd 0.5 1
done




######SCDS warp with deformation constrains, no shape index, init with Id
weight=0.01
mkdir  $scdsWarpFolder/CNT/InitById_$weight

# get the frechet mean CNT contour
txApply -b -i  $cbctFolder/CNT/phase50.mhd -h $scdsWarpFolder/InitById_$weight/scdsWarp_Hfield_0003.mhd -o  $cbctFolder/CNT/atlas_mean_InitById_$weight

for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i  $cbctFolder/CNT/atlas_mean_InitById_$weight.mhd  -h  $scdsWarpFolder/InitById_$weight/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/CNT/InitById_$weight/ref-cbct50-phase$i  
ThresholdImage 3 $scdsWarpFolder/CNT/InitById_$weight/ref-cbct50-phase$i.mhd  $scdsWarpFolder/CNT/InitByI_$weightd/ref-cbct50-phase$i.mhd 0.5 1
done




