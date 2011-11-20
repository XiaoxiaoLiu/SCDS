#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 2
#$ -j y
#$ -pe smp 6



PatientName=$1
PredType=$2
weight=$3
echo $PredType
echo $weight
cbctFolder=/stage/sharonxx/proj/mskcc/$PatientName/cbct
resultsFolder=/stage/sharonxx/proj/mskcc/$PatientName/results

predictionFolder=$resultsFolder/SCDS_predict/$PredType/Hfield
atlasFolder=/stage/sharonxx/proj/mskcc/$PatientName/results/SCDS_warp/$PredType
mkdir /stage/sharonxx/proj/mskcc/$PatientName/results/SCDS_warp/
mkdir $atlasFolder

#initialized by identity transform
if [ $SGE_TASK_ID -eq 2 ]; then
	atlasFolder=$atlasFolder/InitById_$weight
	mkdir $atlasFolder
	~sharonxx/bin/linux/scdsWarp --initWithInputHField=flase  --optWeights  1 $weight $weight $weight    --numberOfImages=6   --inputImageList   $cbctFolder/image/gray/phase*.mhd  --inputHFieldList $predictionFolder/pred*.mhd   --outputImageFilenamePrefix=$atlasFolder/scdsWarp_averageImage_       --outputDeformedImageFilenamePrefix=$atlasFolder/scdsWarp_deformedImage_    --outputHFieldFilenamePrefix=$atlasFolder/scdsWarp_Hfield_   --scaleLevel=4 --numberOfIterations=100  --scaleLevel=2 --numberOfIterations=100   --scaleLevel=1 --numberOfIterations=50  > $atlasFolder/SCDSWarp.log 2>&1
fi


#only use intensity to drive the optimization
if [ $SGE_TASK_ID -eq 3 ]; then
	atlasFolder=$atlasFolder/OnlyIntensity
	mkdir $atlasFolder
	~sharonxx/bin/linux/scdsWarp --initWithInputHField=false  --optWeights  1 0 0 0    --numberOfImages=6   --inputImageList   $cbctFolder/image/gray/phase*.mhd  --inputHFieldList $predictionFolder/pred*.mhd   --outputImageFilenamePrefix=$atlasFolder/scdsWarp_averageImage_       --outputDeformedImageFilenamePrefix=$atlasFolder/scdsWarp_deformedImage_    --outputHFieldFilenamePrefix=$atlasFolder/scdsWarp_Hfield_   --scaleLevel=4 --numberOfIterations=100 --scaleLevel=2 --numberOfIterations=100    --scaleLevel=1 --numberOfIterations=50  > $atlasFolder/SCDSWarp.log 2>&1 
fi
