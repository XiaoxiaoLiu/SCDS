# Xiaoxiao Liu
# scripts for preparing training data for SCDS predction on CBCT
# pat 20 21 25 30, all patients have CBCT and RCCT data, aligned and resampled


# for NCAT, segmentting the CBCT is really easy


# alignment done

# registration using fWarp
registration-fWarp-qsub.bash



# binary seg lung from rcct
## get the lung regions out from the binary threshold

qsub ~sharonxx/scripts/SCDS-cbct/4pats/binarySegRCCT.bash

# filling the holes
qsub ~sharonxx/scripts/SCDS-cbct/4pats/rollingballRCCT.bash







#generate smoothed distance map


#chop the top and bottom slice to make sure closed surface
#matlab/common/makeChoppingMask.m

#for j in 20 21 25 31
#do
#cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/shape
#for i in 00 10 20 30 40 50 60 70 80 90 
#do  
#/home/sharonxx/bin/linux/ImageMath 3  lung-bin-cinePhase$i.mhd  m  ../image/binary-inter/lung-bin-RB-cinePhase$i.mhd  chop_mask.mhd 
#done
#done

#pad image

Y:\scripts\SCDS-cbct\4pats\padBinaryImage.m
#gengerate corresponding particles, run correpondence program

#for j in 20 21 25 31
#do
#cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/shape
#ShapeWorksGroom  groom.params  antialias fastmarching blur
#done


qub ~sharonxx/scripts/SCDS-cbct/4pats/distancemap_qsub.bash



#for i in 00 10 20 30 40 50 60 70 80 90 ; do ImageMath_Ipek  lung-bin-cinePhase$i.DT.mhd   -outfile  lung-bin-cinePhase$i.DT.mhd   -editPixdims 1,1,1 -type float; done
##########################################################
for j in 20 21 25 31
do
ShapeWorksShop shapeWorksRun.params
done


# put seeds on left and right lung
Y:\scripts\SCDS-cbct\4pats\writeSeeds.m


#5.generate corresponding meshes, run postprocessing program
qsub  /home/sharonxx/scripts/SCDS-cbct/4pats/vol2surf.bash
ParticleCorrespondencePostprocessing  --parameterFileName corresMesh.params

# train mesh stats
trainMeshPtsStats_4pats

#segment CBCT
mkdir Z:\proj\mskcc\Patient21\cbct\segmentation
#VC arguments
-p 5 -w 0.1  -n 100 -d 0.1 -i Z:\proj\mskcc\Patient21\cbct\image\gray-inter\phase95.mhd  -s Z:\proj\mskcc\Patient21\rcct\stats\meshPtsStats\meshPtsStats_pat21.txt -m Z:\proj\mskcc\Patient21\rcct\stats\meshPtsStats\meanMesh_pat21.vtk -o  Z:\proj\mskcc\Patient21\cbct\segmentation\scds_phase95.vtk
