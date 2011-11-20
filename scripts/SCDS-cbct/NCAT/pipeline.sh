
#original method, the cbct prediction is not perfectly matched with rcct




run_prediction_inter.m



qsub /home/sharonxx/scripts/SCDS-cbct/NCAT/scdsWarp.bash NCAT3 orig 0.00001


evaluateCNT.bash

#NCAT2 orig
qsub /home/sharonxx/scripts/SCDS-cbct/NCAT/scdsWarp.bash NCAT2 orig 0.00001



########################### old begin ################################
##################################################################

##################################################################
####################### results are not improved:too much differences in NCAT2 and NCAT3, the prediction is too off?? CNT CBCT intensities are good #######################################

#############train on NCAT2 and test NCAT3

#train on NCAT2
NCAT2.bash
evaluateCNT.bash



# test on NCAT3
NCAT3.bash
evaluateCNT.bash


# use NCAT2's SCDS model to test NCAT3 data
		

#segment NCAT3
mkdir /stage/sharonxx/proj/mskcc/NCAT3/cbct/segmentation
mkdir /stage/sharonxx/proj/mskcc/NCAT3/cbct/segmentation/NCAT2_Stats
#
meshDeform-CBCT-NCAT3-NCAT2Stats.bash
#
run_prediction_inter_train_NCAT2_pred_NCAT3.m

qsub /home/sharonxx/scripts/SCDS-cbct/NCAT/scdsWarp.bash NCAT3 NCAT2_stats 0.0001

evaluateCNT.bash

########################### old end ################################
##################################################################




##############randError



#NCAT2
qsub /home/sharonxx/scripts/SCDS-cbct/NCAT/scdsWarp.bash NCAT2 randError_stats 0.0001



evaluateCNT.bash

evaluate_CNT_COG_ref_to_CBCT50.m




#NCAT3
run_prediction_inter_RandError.m

qsub /home/sharonxx/scripts/SCDS-cbct/NCAT/scdsWarp.bash NCAT3 randError_stats 0.0001



evaluateCNT.bash

evaluate_CNT_COG_ref_to_CBCT50.m





########### generate frechet mean of prediction with randomError
PatientName=NCAT2
cbctFolder=/stage/sharonxx/proj/mskcc/$PatientName/cbct
resultsFolder=/stage/sharonxx/proj/mskcc/$PatientName/results
scdsPredictionFolder=$resultsFolder/SCDS_predict/randError_stats
mkdir  $scdsPredictionFolder/deformedImages



for i in 10 20 30 40 50 60 
do txApply -b -i   $cbctFolder/image/gray/phase$i.mhd    -h $scdsPredictionFolder/Hfield/pred-hfield-phase$i.mhd  -o $scdsPredictionFolder/deformedImages/phase$i  
done

~sharonxx/bin/linux/AtlasWerks  --outputImageFilenamePrefix=$scdsPredictionFolder/averageImage_  --outputDeformedImageFilenamePrefix=$scdsPredictionFolder/deformedImage_  --outputHFieldFilenamePrefix=$scdsPredictionFolder/hfield_      --scaleLevel=1  --numberOfIterations=0  $scdsPredictionFolder/deformedImages/phase*.mhd




