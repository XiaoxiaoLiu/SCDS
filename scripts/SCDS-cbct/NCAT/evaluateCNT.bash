

####################################### not working yet################################
##################### CNT quality is bad ############################

# convert CNT to GTV binary image
# align the GTV to the intersected image
itk-snap segmentation

patNo=3
predType=randError_stats;
#NCAT2_stats
weight=0.00001
#orig
# randError_stats;  # 


scriptsFolder=~sharonxx/scripts/SCDS-cbct/patient
dataFolder=/stage/sharonxx/proj/mskcc/NCAT$patNo
cbctFolder=$dataFolder/cbct
cbctImageFolder=$dataFolder/cbct/image

rcctFolder=$dataFolder/rcct
rcctImageFolder=$dataFolder/rcct/image
rcctAtlasFolder=$dataFolder/rcct/atlas
shapeFolder=$dataFolder/rcct/shape
statsFolder=$dataFolder/rcct/stats

resultsFolder=$dataFolder/results

scdsPredictionFolder=$dataFolder/results/SCDS_predict/$predType
scdsWarpFolder=$dataFolder/results/SCDS_warp/$predType








mkdir $scdsPredictionFolder/CNT
mkdir  $scdsWarpFolder/CNT

######SCDS warp with deformation constrains, no shape index, init with Id
mkdir  $scdsWarpFolder/CNT/InitById_$weight

# get the frechet mean CNT contour
txApply -b -i  $cbctFolder/CNT/phase40.mhd -h $scdsWarpFolder/InitById_$weight/scdsWarp_Hfield_0003.mhd -o  $cbctFolder/CNT/atlas_mean_InitById_$weight

for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i  $cbctFolder/CNT/atlas_mean_InitById_$weight.mhd  -h  $scdsWarpFolder/InitById_$weight/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/CNT/InitById_$weight/ref-cbct50-phase$i  
ThresholdImage 3 $scdsWarpFolder/CNT/InitById_$weight/ref-cbct50-phase$i.mhd  $scdsWarpFolder/CNT/InitById_$weight/ref-cbct50-phase$i.mhd 0.5 1
done















##### SCDS prediction


# get the frechet mean CNT contour
txApply -b -i  $cbctFolder/CNT/phase40.mhd -h $scdsPredictionFolder/Hfield/pred-hfield-phase40.mhd -o  $scdsPredictionFolder/CNT/atlas_mean_pred

for i in 10 20 30 40 50 60 
do txApply -f -i   $scdsPredictionFolder/CNT/atlas_mean_pred.mhd  -h $scdsPredictionFolder/Hfield/pred-hfield-phase$i.mhd  -o $scdsPredictionFolder/CNT/ref-cbct50-phase$i  
ThresholdImage 3  $scdsPredictionFolder/CNT/ref-cbct50-phase$i.mhd $scdsPredictionFolder/CNT/ref-cbct50-phase$i.mhd 0.5 1
done











#####SCDS warp OnlyIntensity
mkdir  $scdsWarpFolder/CNT/OnlyIntensity

# get the frechet mean CNT contour
txApply -b -i  $cbctFolder/CNT/phase40.mhd -h $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_0003.mhd -o  $cbctFolder/CNT/atlas_mean_Int

for i in 0000 0001 0002 0003 0004 0005 
do txApply -f -i    $cbctFolder/CNT/atlas_mean_Int.mhd  -h  $scdsWarpFolder/OnlyIntensity/scdsWarp_Hfield_$i.mhd  -o $scdsWarpFolder/CNT/OnlyIntensity/ref-cbct50-phase$i  
ThresholdImage 3 $scdsWarpFolder/CNT/OnlyIntensity/ref-cbct50-phase$i.mhd $scdsWarpFolder/CNT/OnlyIntensity/ref-cbct50-phase$i.mhd 0.5 1
done






# produce the  CNT using interpolated Hfield, copy RCCT CNT cinePhase20.mhd to replace CBCT CNT phase20.mhd 
txApply -b -i  $cbctFolder/CNT/phase20.mhd -h $cbctFolder/interpolatedHfield/hfield_2.mhd -o  $cbctFolder/CNT/atlas_mean_interpH
mkdir  $cbctFolder/interpolatedHfield/CNT
for i in 1 2 3 4 5 6 
do txApply -f -i  $cbctFolder/CNT/atlas_mean_interpH.mhd  -h  $cbctFolder/interpolatedHfield/hfield_$i.mhd  -o  $cbctFolder/interpolatedHfield/CNT/ref-cbct50-phase$i  
ThresholdImage 3 $cbctFolder/interpolatedHfield/CNT/ref-cbct50-phase$i.mhd  $cbctFolder/interpolatedHfield/CNT/ref-cbct50-phase$i.mhd 0.5 1
done



