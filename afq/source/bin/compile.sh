#!/bin/bash

cat > build.m <<END

addpath(genpath('/data/localhome/glerma/soft/RTP-pipeline'));
rmpath(genpath('/data/localhome/glerma/soft/RTP-pipeline/local'));

addpath(genpath('/data/localhome/glerma/toolboxes/jsonlab'));
addpath(genpath('/data/localhome/glerma/toolboxes/JSONio'));

addpath(genpath('/data/localhome/glerma/soft/encode'));
addpath(genpath('/data/localhome/glerma/soft/app-life'));


addpath(genpath('/data/localhome/glerma/toolboxes/freesurfer_mrtrix_afni_matlab_tools'));


mcc -m -R -nodisplay -a /data/localhome/glerma/soft/RTP-pipeline/afq/includeFiles -a /data/localhome/glerma/soft/encode/mexfiles  -d compiled RTP.m

exit
END
/software/matlab/r2018b/bin/matlab -nodisplay -nosplash -r build && rm build.m


# The compiled file is bigger than 100Mb, then it fails when pushing to github
# Use the command below to delete it from all history
# Add a gitignore so that it never goes back anything in the /compiled folder
# git filter-branch --force --index-filter 'git rm --cached --ignore-unmatch /data/localhome/glerma/soft/afq-pipeline/afq/source/bin/compiled/AFQ_StandAlone_QMR' --prune-empty --tag-name-filter cat -- --all









