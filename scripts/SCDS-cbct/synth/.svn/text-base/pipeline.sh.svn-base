
#genearte synthetic data

matlab/sythData/spheres4DGen.m




#use frechet mean image as the reference image

#training
cd /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1

mkdir ./atlas
AtlasWerks --outputImageFilenamePrefix=./atlas/averageImage_  --outputDeformedImageFilenamePrefix=./atlas/deformedImage_  --outputHFieldFilenamePrefix=./atlas/hfield_  --scaleLevel=1  --numberOfIterations=500  ./orig/sphere*.mhd   > ./atlas/SCDS.log 2>&1


#run statsitsics
run_pred_inter.m



# run test1
# 1) sanity check
scdsWarp.bash
optimizationCurve.m

# 2) random effects adds to the correlation
scdsWarp_randCorr.bash
optimizationCurve_randCorr.m





























##############################################  no longer used #############################################
######################### test2 #####################

cd /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test2 
 
## ground truth
mkdir ./atlas
AtlasWerks --outputImageFilenamePrefix=./atlas/averageImage_  --outputDeformedImageFilenamePrefix=./atlas/deformedImage_  --outputHFieldFilenamePrefix=./atlas/hfield_  --scaleLevel=1  --numberOfIterations=500  ./orig/sphere*.mhd   > ./atlas/SCDS.log 2>&1






####comparison 
cd /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test2    
mkdir atlas


#single scale only with intensity
mkdir atlas/singleScale_Int/

scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  100 0 0 0    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ./pred/pred-hfield*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./atlas/singleScale_Int/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale_Int/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale_Int/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale_Int/SCDS.log 2>&1

mkdir atlas/singleScale_H/
scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  0 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ./pred/pred-hfield*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./atlas/singleScale_H/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale_H/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale_H/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale_H/SCDS.log 2>&1

#single scale 
mkdir atlas/singleScale/
scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  100 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd  --inputHFieldList ./pred/pred-hfield*.mhd   --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./atlas/singleScale/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale/SCDS.log 2>&1
  
 

## single scale with shape score
mkdir atlas/singleScale_S

# the shape score coming from the stats_P.score(:,1): the first PC score
scdsWarp  --writeZSlice=32  --initWithInputHField=false --optWeights  100 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ./pred/pred-hfield*.mhd    --inputCompHFieldList ./atlas/hfield_*.mhd   --shapeScores  0.9400    0.5587    0.0322    0.5818    1.0000    1.0000    0.5818    0.0322    0.5587    0.9400 --outputImageFilenamePrefix=./atlas/singleScale_S/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale_S/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale_S/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale_S/SCDS.log 2>&1




######################### test3 #####################

cd /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test3
 
## ground truth
mkdir ./atlas
AtlasWerks --outputImageFilenamePrefix=./atlas/averageImage_  --outputDeformedImageFilenamePrefix=./atlas/deformedImage_  --outputHFieldFilenamePrefix=./atlas/hfield_  --scaleLevel=1  --numberOfIterations=500  ./orig/sphere*.mhd   > ./atlas/SCDS.log 2>&1



#run statsitsics
run_pred_inter.m







####comparison 
cd /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1    
mkdir atlas


#single scale only with intensity
mkdir atlas/singleScale_Int/

scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  100 0 0 0    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ./pred/pred-hfield*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./atlas/singleScale_Int/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale_Int/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale_Int/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale_Int/SCDS.log 2>&1

mkdir atlas/singleScale_H/
scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  0 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ../test1/atlas/hfield*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./atlas/singleScale_H/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale_H/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale_H/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale_H/SCDS.log 2>&1

#single scale 
mkdir atlas/singleScale/
scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  100 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd  --inputHFieldList  ../test1/atlas/hfield*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./atlas/singleScale/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale/SCDS.log 2>&1
  
 

## single scale with shape score
mkdir atlas/singleScale_S

# the shape score coming from the stats_P.score(:,1): the first PC score
scdsWarp  --writeZSlice=32  --initWithInputHField=false --optWeights  100 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ./pred/pred-hfield*.mhd    --inputCompHFieldList ./atlas/hfield_*.mhd   --shapeScores  1.0000    0.6181         0    0.6181    1.0000    1.0000    0.6181         0    0.6181    1.0000 --outputImageFilenamePrefix=./atlas/singleScale_S/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./atlas/singleScale_S/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./atlas/singleScale_S/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./atlas/singleScale_S/SCDS.log 2>&1


