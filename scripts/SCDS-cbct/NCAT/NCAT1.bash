# for NCAT data, no prealignment needed




NUM=1;
spacing=(1.52 1.52 1.52)  # CBCT's image sapcing for resampling purpose

#default data storgage settings
scriptsFolder=~sharonxx/scripts/SCDS-cbct/NCAT
dataFolder=/stage/sharonxx/proj/mskcc/NCAT$NUM
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

mkdir $shapeFolder
mkdir $statsFolder
mkdir $resultsFolder
mkdir $scdsPredictionFolder
mkdir $scdsWarpFolder





#resample rcct to match the cbct
cd /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/image
mkdir gray-resample
for i in 00 10 20 30 40 50 60 70 80 90 
do ResampleImageBySpacing 3 gray/cinePhase$i.mhd gray-resample/cinePhase$i.mhd 0.742 0.742 1.52
done

#matlab function /common
intersectionRegionExtract_NCAT.m



#2. run registration  
#registration-fWarp-qsub-ncat.bash

qsub run_atlasWerks.bash



#3.segmentation of the rcct

NUM=3;  #with tumor
cd /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/image
mkdir  /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/image/binary-inter
for i in 00 10 20 30 40 50 60 70 80 90 
do 
ThresholdImage 3 ./gray-inter/cinePhase$i.mhd ./binary-inter/lung-bin-cinePhase$i.mhd  550 1800
RollingBall  ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd  0 0  2  2  2  2  50
ConnectedThreshold ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd   1  2  139 231 45  401 231 34 
done





#4. Build the PDMS with correpondence

cd /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/shape


#################### tbd ##################
# ???something new for NCAT2!!!!
#resample the image into [2 2 2]
##for i in 00 10 20 30 40 50 60 70 80 90;do  ResampleImageBySpacing  3 lung-bin-cinePhase$i.mhd  lung-bin-cinePhase$i.mhd  1 1 1; done

##note that if the images are too big, the particle system will be very slow
#cd /stage/sharonxx/proj/mskcc/NCAT2/rcct/image/binary-inter
#for i in 00 10 20 30 40 50 60 70 80 90;do  ResampleImageBySpacing  3 lung-bin-cinePhase$i.mhd  ../../shape/lung-bin-cinePhase$i.mhd  2 2 2; done


##pad image, put binary image in the shape folder, 
#scripts/NCAT/padBinaryImage.m
################  tbd ###################


# Crop the binary image and pad the boundaires
crop=(80 110 0  370 280 100)   #(starting index,  size)parameter for CropAndPad
for i in  00 10 20 30 40 50 60 70 80 90 ; 
do CropAndPad  $rcctImageFolder/binary-inter/lung-bin-cinePhase$i.mhd  $shapeFolder/lung-bin-crop-pad-cinePhase$i.mhd  ${crop[0]}  ${crop[1]}  ${crop[2]} ${crop[3]} ${crop[4]} ${crop[5]}  8 8 8 ;
done



qsub  /home/sharonxx/scripts/SCDS-cbct/NCAT/distancemap_qsub-ncat.bash
qsub  /home/sharonxx/scripts/SCDS-cbct/NCAT/vol2surf-nact.bash


./writeSeeds_NCAT.m 


Correspondecne on XYZ: {
#tuning the correspondence weight around 10
mkdir xyz
cd xyz
Correspondence particle.params
#5. train SCDS statistics
trainPtsStats.m
}


~ipek/NeuroLib-build/bin/ParticleCorrespondencePostprocessingTPS --parameterFileName corresMesh.params
#################################################################################################
trainMeshPtsStats.m



 qsub ~sharonxx/scripts/SCDS-cbct/NCAT/meshDeform-CBCT.bash


# evaluate the fitting results , comparing with  RCCT  shape PC scores
evalDeformResults.m


#8. predict CBCT displacement vector fields
run_predict_inter.m





#9.producing average image
NUM=1;
mkdir /stage/sharonxx/proj/mskcc/NCAT$NUM/results/predict-CCA_XY/deformedImage
for n in 10 20 30 40 50 60 
do
txApply -b -i /stage/sharonxx/proj/mskcc/NCAT$NUM/cbct/image/gray/phase$n.mhd -h /stage/sharonxx/proj/mskcc/NCAT$NUM/results/predict-CCA_XY/Hfield-ALL/pred-hfield-phase$n.mhd  -o /stage/sharonxx/proj/mskcc/NCAT$NUM/results/predict-CCA_XY/deformedImage/deformedPhase$n.mhd
done

#diapghram comparison
for n in 10 20 30 40 50 60 
do
txApply -b -i /stage/sharonxx/proj/mskcc/NCAT$NUM/cbct/image/gray/phase$n.mhd -h /stage/sharonxx/proj/mskcc/NCAT$NUM/results/predict-MLR_XY/Hfield-ALL/dg-pred-hfield-phase$n.mhd  -o /stage/sharonxx/proj/mskcc/NCAT$NUM/results/predict-MLR_XY/deformedImage/dg-deformedPhase$n.mhd
done
#generate the mean image
scdsWarp  --outputImageFilenamePrefix=   --outputDeformedImageFilenamePrefix=   --outputHFieldFilenamePrefix=

Image
