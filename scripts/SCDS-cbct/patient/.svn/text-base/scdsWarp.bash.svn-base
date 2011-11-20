#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 2
#$ -j y
#$ -pe smp 6



PatientName=$1;
weight=$2

cbctFolder=/stage/sharonxx/proj/mskcc/$PatientName/cbct
resultsFolder=/stage/sharonxx/proj/mskcc/$PatientName/results
atlasFolder=/stage/sharonxx/proj/mskcc/$PatientName/results/SCDS_warp
mkdir $atlasFolder

#initialized by hfield
if [ $SGE_TASK_ID -eq 1 ]; then
	atlasFolder=/stage/sharonxx/proj/mskcc/$PatientName/results/SCDS_warp/InitByH_$weight
	mkdir $atlasFolder
	~sharonxx/bin/linux/scdsWarp --initWithInputHField=true  --optWeights  1 $weight  $weight $weight    --numberOfImages=6   --inputImageList   $cbctFolder/image/gray-inter/phase*.mhd  --inputHFieldList $resultsFolder/SCDS_predict/Hfield/pred*.mhd   --outputImageFilenamePrefix=$atlasFolder/scdsWarp_averageImage_       --outputDeformedImageFilenamePrefix=$atlasFolder/scdsWarp_deformedImage_    --outputHFieldFilenamePrefix=$atlasFolder/scdsWarp_Hfield_   --scaleLevel=4 --numberOfIterations=100 --scaleLevel=2 --numberOfIterations=100   --scaleLevel=1 --numberOfIterations=50 > $atlasFolder/SCDSWarp.log 2>&1
fi

#initialized by identity transform
if [ $SGE_TASK_ID -eq 2 ]; then
	atlasFolder=/stage/sharonxx/proj/mskcc/$PatientName/results/SCDS_warp/InitById_$weight
	mkdir $atlasFolder
	~sharonxx/bin/linux/scdsWarp --initWithInputHField=flase  --optWeights  1 $weight  $weight $weight   --numberOfImages=6   --inputImageList   $cbctFolder/image/gray-inter/phase*.mhd  --inputHFieldList $resultsFolder/SCDS_predict/Hfield/pred*.mhd   --outputImageFilenamePrefix=$atlasFolder/scdsWarp_averageImage_       --outputDeformedImageFilenamePrefix=$atlasFolder/scdsWarp_deformedImage_    --outputHFieldFilenamePrefix=$atlasFolder/scdsWarp_Hfield_   --scaleLevel=4 --numberOfIterations=100 --scaleLevel=2 --numberOfIterations=100   --scaleLevel=1 --numberOfIterations=50 > $atlasFolder/SCDSWarp.log 2>&1
fi


#only use intensity to drive the optimization
if [ $SGE_TASK_ID -eq 3 ]; then
	atlasFolder=/stage/sharonxx/proj/mskcc/$PatientName/results/SCDS_warp/OnlyIntensity
	mkdir $atlasFolder
	~sharonxx/bin/linux/scdsWarp --initWithInputHField=false  --optWeights  1 0 0 0    --numberOfImages=6   --inputImageList   $cbctFolder/image/gray-inter/phase*.mhd  --inputHFieldList $resultsFolder/SCDS_predict/Hfield/pred*.mhd   --outputImageFilenamePrefix=$atlasFolder/scdsWarp_averageImage_       --outputDeformedImageFilenamePrefix=$atlasFolder/scdsWarp_deformedImage_    --outputHFieldFilenamePrefix=$atlasFolder/scdsWarp_Hfield_   --scaleLevel=4 --numberOfIterations=100 --scaleLevel=2 --numberOfIterations=100   --scaleLevel=1 --numberOfIterations=25 > $atlasFolder/SCDSWarp.log  2>&1
fi
