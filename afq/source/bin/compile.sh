#!/bin/bash
# module load matlab/2017a

cat > build.m <<END

addpath(genpath('~/soft/AFQ'));
addpath(genpath('~/soft/RTP-pipeline'));
addpath(genpath('~/soft/vistasoft'));
addpath(genpath('~/soft/spm8'));
rmpath(genpath('~/soft/spm8/toolbox/Beamforming'));
rmpath(genpath('~/soft/spm8/toolbox/DARTEL'));
rmpath(genpath('~/soft/spm8/toolbox/MEEGtools'));
rmpath(genpath('~/soft/spm8/toolbox/FieldMap'));
rmpath(genpath('~/soft/spm8/toolbox/dcm_meeg'));
rmpath(genpath('~/soft/spm8/external/bemcp'));
rmpath(genpath('~/soft/spm8/external/ctf'));
rmpath(genpath('~/soft/spm8/external/eeprobe'));
rmpath(genpath('~/soft/spm8/external/fieldtrip'));
rmpath(genpath('~/soft/spm8/external/mne'));
rmpath(genpath('~/soft/spm8/external/yokogawa'));

addpath(genpath('~/soft/jsonlab'));
addpath(genpath('~/soft/encode'));
addpath(genpath('~/soft/JSONio'));
addpath(genpath('~/soft/app-life'));
addpath(genpath('~/soft/freesurfer_mrtrix_afni_matlab_tools'));


mcc -m -R -nodisplay -a /data/localhome/glerma/soft/RTP-pipeline/afq/includeFiles -a /data/localhome/glerma/soft/encode/mexfiles -a /data/localhome/glerma/soft/vistasoft/mrDiffusion -a /data/localhome/glerma/soft/AFQ/templates/labelMaps -d compiled AFQ_StandAlone_QMR.m

exit
END
Matlabr2017a -nodisplay -nosplash -r build && rm build.m


# Teh compiled file is bigger than 100Mb, then it fails when pushing to github
# Use the command below to delete it from all history
# Add a gitignore so that it never goes back anything in the /compiled folder
# git filter-branch --force --index-filter 'git rm --cached --ignore-unmatch /data/localhome/glerma/soft/afq-pipeline/afq/source/bin/compiled/AFQ_StandAlone_QMR' --prune-empty --tag-name-filter cat -- --all

