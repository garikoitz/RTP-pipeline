#!/bin/bash

while getopts "b:d:p:t:h:u:g:m:q:c:a:r:x:f:" opt;do
	case $opt in
		b) BASEDIR="$OPTARG";;
		d) codedir="$OPTARG";;
		p) MNI_package="$OPTARG";;
		t) MNI_template="$OPTARG";;
		h) host="$OPTARG";;
		u) qsub="$OPTARG";;
		g) logdir="$OPTARG";;
		m) mem="$OPTARG";;
		q) que="$OPTARG";;
		c) core="$OPTARG";;
		a) ants_ver="$OPTARG";;
		r) mrtrix_ver="$OPTARG";;
		x) gcc_ver="$OPTARG";;
		f) fsl_ver="$OPTARG";;
	esac
done	

module load $gcc_ver
module load $mrtrix_ver
module load $fsl_ver

# for loop for subjects in the DB
for proj in $(cat ${BASEDIR}/subSesList.txt);do

printf "Working on $proj \n"
MRI_DIR=${BASEDIR}/${proj}/structural

printf "Performing skull stripping \n"
mrtransform -force -template ${BASEDIR}/${proj}/structural/t1w_acpc.nii.gz ${BASEDIR}/${proj}/dt6/bin/brainMask.nii.gz  ${BASEDIR}/${proj}/structural/brainmask_out.nii.gz

fslmaths ${BASEDIR}/${proj}/structural/t1w_acpc.nii.gz -mul ${BASEDIR}/${proj}/structural/brainmask_out.nii.gz ${BASEDIR}/${proj}/structural/brain.nii.gz

done
