#$ -S /bin/bash
#$ -o /home/sharonxx/tmp/
#$ -t 1
#$ -j y
# registration using fWarp


NUM=2;
echo "this is job  ncat phase$i"

mkdir /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct/atlas
cd  /stage/sharonxx/proj/mskcc/NCAT$NUM/rcct

/home/sharonxx/bin/linux/AtlasWerks --outputImageFilenamePrefix=./atlas/averageImage_  --outputDeformedImageFilenamePrefix=./atlas/deformedImage_  --outputHFieldFilenamePrefix=./atlas/hfield_   --scaleLevel=4  --numberOfIterations=100   --scaleLevel=2  --numberOfIterations=100 --scaleLevel=1  --numberOfIterations=50  ./image/gray-inter/cinePhase*.mhd   > ./atlas/atlasWerks.log 2>&1
