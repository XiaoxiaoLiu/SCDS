#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1-1
#$ -j y
#$ -pe smp 10

patName=$1
############## Training atlas 
cd /stage/sharonxx/proj/mskcc/$patName/rcct
mkdir ./atlas
/home/sharonxx/bin/linux/AtlasWerks --outputImageFilenamePrefix=./atlas/averageImage_  --outputDeformedImageFilenamePrefix=./atlas/deformedImage_  --outputHFieldFilenamePrefix=./atlas/hfield_   --scaleLevel=4  --numberOfIterations=200  --scaleLevel=2  --numberOfIterations=200  --scaleLevel=1  --numberOfIterations=100  ./image/gray-inter/cinePhase*0.mhd    > ./atlas/atlasWerks.log 2>&1


#############  Test Intensity atlas on CBCT

cd /stage/sharonxx/proj/mskcc/$patName/cbct
mkdir ./atlas
/home/sharonxx/bin/linux/AtlasWerks --outputImageFilenamePrefix=./atlas/averageImage_  --outputDeformedImageFilenamePrefix=./atlas/deformedImage_  --outputHFieldFilenamePrefix=./atlas/hfield_    --scaleLevel=4  --numberOfIterations=200  --scaleLevel=2  --numberOfIterations=200  --scaleLevel=1  --numberOfIterations=100   ./image/gray-inter/phase*.mhd    >  ./atlas/atlasWerks.log 2>&1 
