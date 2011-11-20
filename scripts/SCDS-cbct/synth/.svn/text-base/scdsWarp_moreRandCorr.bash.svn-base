#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-3
#$ -j y


cd /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1    

postfix=randCorr-more
atlasFolder=atlas-$postfix 

mkdir $atlasFolder


if [ $SGE_TASK_ID = "1" ]; then
#only Intensity
mkdir $atlasFolder/singleScale_Int
~sharonxx/bin/linux/scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  1 0 0 0    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ./pred-$postfix/pred-hfield*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./$atlasFolder/singleScale_Int/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./$atlasFolder/singleScale_Int/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./$atlasFolder/singleScale_Int/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./$atlasFolder/singleScale_Int/SCDS.log 2>&1
fi




if [ $SGE_TASK_ID = "2" ]; then

#only use H
mkdir $atlasFolder/singleScale_H/
~sharonxx/bin/linux/scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  0 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList  ./pred-$postfix/pred-hfield*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./$atlasFolder/singleScale_H/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./$atlasFolder/singleScale_H/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./$atlasFolder/singleScale_H/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./$atlasFolder/singleScale_H/SCDS.log 2>&1
fi



if [ $SGE_TASK_ID = "3" ]; then
#combine H and int
mkdir $atlasFolder/singleScale/
~sharonxx/bin/linux/scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  100 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd  --inputHFieldList  ./pred-$postfix/pred-hfield*.mhd   --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./$atlasFolder/singleScale/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./$atlasFolder/singleScale/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./$atlasFolder/singleScale/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./$atlasFolder/singleScale/SCDS.log 2>&1
 fi





