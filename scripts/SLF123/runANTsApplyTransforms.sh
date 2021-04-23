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

for ROI in $roidir/*.nii.gz;do
roiname=$(basename $ROI)
roiname=${roiname%%.*}

printf "Transforming $roiname \n"

antsApplyTransforms -d 3 \
-i $ROI \
-r $ref \
-n BSpline \
-t $transform1 \
-t $transform2 \
-o $roidir/${roiname}-MNISegment.nii.gz

done
