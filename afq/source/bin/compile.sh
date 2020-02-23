#!/bin/bash

cat > build.m <<END

% addpath(genpath('/data/localhome/lmx/TESTDATA/AFQ'));
% addpath(genpath('/data/localhome/lmx/TESTDATA/RTP-pipeline'));
addpath(genpath('/data/localhome/glerma/soft/rtppipelinemlx/RTP-pipeline'));
addpath(genpath('/data/localhome/lmx/TESTDATA/vistasoft'));
% addpath(genpath('/data/localhome/glerma/soft/spm8'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/toolbox/Beamforming'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/toolbox/DARTEL'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/toolbox/MEEGtools'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/toolbox/FieldMap'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/toolbox/dcm_meeg'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/external/bemcp'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/external/ctf'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/external/eeprobe'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/external/fieldtrip'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/external/mne'));
% rmpath(genpath('/data/localhome/glerma/soft/spm8/external/yokogawa'));

addpath(genpath('/data/localhome/glerma/soft/jsonlab'));
addpath(genpath('/data/localhome/glerma/soft/encode'));
addpath(genpath('/data/localhome/glerma/soft/JSONio'));
addpath(genpath('/data/localhome/glerma/soft/app-life'));
addpath(genpath('/data/localhome/glerma/soft/freesurfer_mrtrix_afni_matlab_tools'));


mcc -m -R -nodisplay -a /data/localhome/lmx/TESTDATA/RTP-pipeline/afq/includeFiles -a /data/localhome/glerma/soft/encode/mexfiles -a /data/localhome/lmx/TESTDATA/vistasoft/mrDiffusion -a /data/localhome/lmx/TESTDATA/RTP-pipeline/afq/repo/templates/labelMaps -d compiled RTP.m

exit
END
/software/matlab/r2018b/bin/matlab -nodisplay -nosplash -r build && rm build.m


# Teh compiled file is bigger than 100Mb, then it fails when pushing to github
# Use the command below to delete it from all history
# Add a gitignore so that it never goes back anything in the /compiled folder
# git filter-branch --force --index-filter 'git rm --cached --ignore-unmatch /data/localhome/glerma/soft/afq-pipeline/afq/source/bin/compiled/AFQ_StandAlone_QMR' --prune-empty --tag-name-filter cat -- --all

