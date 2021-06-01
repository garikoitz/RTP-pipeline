#!/bin/bash

while getopts "b:f:" opt;do
	case $opt in
		b) BASEDIR="$OPTARG";;	
		f) fsl_ver="$OPTARG";;
	esac
done	

module load $fsl_ver

# for loop for subjects in the DB
for proj in $(cat ${BASEDIR}/subjectList.txt);do

printf "Working on $proj \n"

for ROI in $(cat ${BASEDIR}/roisList.txt);do

printf "Binarizing $ROI $"
fslmaths $BASEDIR/$proj/ROIs/$ROI-MNISegment-linear.nii.gz -bin $BASEDIR/$proj/ROIs/$ROI-MNISegment-linear-binary.nii.gz

done
done
