


#training
#1. build atlas, generating hfields   ----OnlyWriteFinalResults=true  not working!!!!!!!!!
mkdir /stage/sharonxx/proj/mskcc/NCAT1/rcct/atlas
AtlasWerks  --OnlyWriteFinalResults=true --outputImageFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT1/rcct/atlas/averageImage_  --outputDeformedImageFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT1/rcct/atlas/deformedImage_  --outputHFieldFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT1/rcct/atlas/hfield_  --scaleLevel=4  --numberOfIterations=100  --scaleLevel=2  --numberOfIterations=50  --scaleLevel=1  --numberOfIterations=25  /stage/sharonxx/proj/mskcc/NCAT1/rcct/image/gray-inter/cinePhase*0.mhd    > /stage/sharonxx/proj/mskcc/NCAT1/rcct/atlas/atlasWerks.log 2>&1

#2 training statistics





#testing

#sanity check, need to change the inputHfields
mkdir /stage/sharonxx/proj/mskcc/NCAT1/cbct_atlas
scdsWarp  --writeZSlice=50 --initWithInputHField=true --optWeights  1 1 1 1     --numberOfImages=6   --inputImageList   /stage/sharonxx/proj/mskcc/NCAT1/cbct/image/gray/phase*.mhd      --inputHFieldList /stage/sharonxx/proj/mskcc/NCAT1/prediction/predict-lung/Hfield-ALL/pred-hfield*.mhd    --outputImageFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT1/cbct_atlas/averageImage_       --outputDeformedImageFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT1/cbct_atlas/deformedImage_    --outputHFieldFilenamePrefix=/stage/sharonxx/proj/mskcc/NCAT1/cbct_atlas/deformationField_    --scaleLevel=4 --numberOfIterations=25     > /stage/sharonxx/proj/mskcc/NCAT1/cbct_atlas/SCDS.log 2>&1




Z:\proj\mskcc\NCAT1\prediction\predict-lung\Hfield-ALL


#vc debug

working direc:

--writeZSlice=50  --optWeights  1 1 1 1    --numberOfImages=2   --inputImageList   z:/proj/mskcc/NCAT1/cbct/image/gray/phase10.mhd   z:/proj/mskcc/NCAT1/cbct/image/gray/phase20.mhd    --inputHFieldList z:/proj/mskcc/NCAT1/prediction/predict-lung/Hfield-ALL/pred-hfield-phase10.mhd z:/proj/mskcc/NCAT1/prediction/predict-lung/Hfield-ALL/pred-hfield-phase20.mhd   --outputImageFilenamePrefix=z:/proj/mskcc/NCAT1/cbct/atlas/scdsWarp/averageImage_       --outputDeformedImageFilenamePrefix=z:/proj/mskcc/NCAT1/cbct/atlas/scdsWarp/deformedImage_    --outputHFieldFilenamePrefix=z:/proj/mskcc/NCAT1/cbct/atlas/scdsWarp/deformationField_     --scaleLevel=1 --numberOfIterations=25

