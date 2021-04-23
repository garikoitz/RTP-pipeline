#!/bin/bash

# setting variables
while getopts "b:d:p:t:h:u:g:m:q:c:a:" opt;do
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
	esac
done	

# look for the MNI template folder, if it does not exist, we create the folder and unzip the template
if [[ ! -d ${BASEDIR}/$MNI_package ]]; then
if [[ -f ${BASEDIR}/${MNI_package}.zip ]]; then
mkdir ${BASEDIR}/$MNI_package
unzip -d ${BASEDIR}/$MNI_package ${BASEDIR}/${MNI_package}.zip
else
printf "Please, provide the ${MNI_package}.zip template \n" 
exit 1
fi
fi

# for loop for subjects in the DB
for proj in $(cat ${BASEDIR}/subjectList.txt);do

MRI_DIR=${BASEDIR}/${proj}/structural
ROI_DIR=${BASEDIR}/${proj}/ROIs

if [ "$qsub" == "true" ];then

if [ "$host" == "BCBL" ];then	
qsub -q $que -N sub_${proj}_transform \
     -o ${logdir}/${proj}_transform.o \
     -e ${logdir}/${proj}_transform.e \
     -l mem_free=$mem \
     -v roidir=${ROI_DIR},MRI_DIR=${MRI_DIR},ants_ver=$ants_ver \
     ${codedir}/runANTsApplyTransforms.sh	
fi

if [ "$host" == "DIPC" ];then

qsub -q $que -l mem=$mem,nodes=1:ppn=$core \
	-N sub_$[proj]_transform \
	-o ${logdir}/${proj}_transform.o \
	-e ${logdir}/${proj}_transform.e \
	-v roidir=${ROI_DIR},MRI_DIR=${MRI_DIR},ants_ver=$ants_ver \
	${codedir}/runANTsApplyTransforms.sh
fi

else
antsApplyTransforms -d 3 \
-i ${ROI} \
-r ${MRI_DIR}/antsWarped.nii.gz \
-n BSpline \
-t ${MRI_DIR}/ants1Warp.nii.gz \
-t ${MRI_DIR}/ants0GenericAffine.mat \
-o ${ROI_DIR}/${roiname}-MNISegment.nii.gz
fi

done

