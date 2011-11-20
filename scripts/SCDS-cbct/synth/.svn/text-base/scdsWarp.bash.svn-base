#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 4
#$ -j y


cd /stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1    


atlasFolder=atlas

mkdir $atlasFolder


if [ $SGE_TASK_ID = "1" ]; then
#only Intensity
mkdir $atlasFolder/singleScale_Int
~sharonxx/bin/linux/scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  1 0 0 0    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList ./atlas/hfield_*.mhd   --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./$atlasFolder/singleScale_Int/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./$atlasFolder/singleScale_Int/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./$atlasFolder/singleScale_Int/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./$atlasFolder/singleScale_Int/SCDS.log 2>&1
fi




if [ $SGE_TASK_ID = "2" ]; then

#only use H
mkdir $atlasFolder/singleScale_H/
~sharonxx/bin/linux/scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  0 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd    --inputHFieldList  ./atlas/hfield_*.mhd  --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./$atlasFolder/singleScale_H/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./$atlasFolder/singleScale_H/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./$atlasFolder/singleScale_H/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./$atlasFolder/singleScale_H/SCDS.log 2>&1
fi



if [ $SGE_TASK_ID = "3" ]; then
#combine H and int
mkdir $atlasFolder/singleScale/
~sharonxx/bin/linux/scdsWarp  --writeZSlice=32 --initWithInputHField=false --optWeights  100 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd  --inputHFieldList ./atlas/hfield_*.mhd   --inputCompHFieldList ./atlas/hfield_*.mhd  --outputImageFilenamePrefix=./$atlasFolder/singleScale/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./$atlasFolder/singleScale/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./$atlasFolder/singleScale/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./$atlasFolder/singleScale/SCDS.log 2>&1
 fi



if [ $SGE_TASK_ID = "4" ]; then
#with shape_Score
mkdir $atlasFolder/singleScale_S
# the shape score coming from the stats_P.score(:,1): the first PC score
~sharonxx/bin/linux/scdsWarp  --writeZSlice=32  --initWithInputHField=false --optWeights  100 1 1 1    --numberOfImages=10 --inputImageList  ./noisy/sphere_noisy1_*.mhd     --inputHFieldList ./atlas/hfield_*.mhd   --inputCompHFieldList ./atlas/hfield_*.mhd   --shapeScores 0.9180    0.6155    0.0494    0.5829    1.0000    1.0000    0.5829    0.0494    0.6155 0.9180  --outputImageFilenamePrefix=./$atlasFolder/singleScale_S/SCDS_averageImage_  --outputDeformedImageFilenamePrefix=./$atlasFolder/singleScale_S/SCDS_deformedImage_  --outputHFieldFilenamePrefix=./$atlasFolder/singleScale_S/SCDS_hfield_  --scaleLevel=1 --numberOfIterations=300  > ./$atlasFolder/singleScale_S/SCDS.log 2>&1
fi

