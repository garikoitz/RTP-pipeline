#!/bin/bash
module load $ants_ver

date;

ref=${MRI_DIR}/antsWarped.nii.gz
transform1=${MRI_DIR}/ants1Warp.nii.gz
transform2=${MRI_DIR}/ants0GenericAffine.mat

printf "#### roi directory: $roidir \n"
printf "#### reference image: $ref \n"
printf "#### transform1: $transform1 \n"
printf "#### transform2: $transform2 \n"
printf "#### BASEDIR: $BASEDIR \n"

for ROI in $(cat $BASEDIR/roisList.txt);do

printf "Transforming $ROI \n"

antsApplyTransforms -d 3 \
-i $roidir/$ROI.nii.gz \
-r $ref \
-n linear \
-t $transform1 \
-t $transform2 \
-o $roidir/$ROI-MNISegment-linear.nii.gz

done
