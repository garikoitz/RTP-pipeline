%% Initialize 
% Set up vistasoft and probabilisticMaps functions and your basedir
% -------------------------------------------------------------------------
% load vistasoft functions
addpath(genpath('/export/home/llecca/llecca/Scripts/SLF123/vistasoft'));
% load probabilisticMap functions
addpath '/export/home/llecca/public/Gari/SLF123/codes/probabilisticMaps'
% define subjects database
BASEDIR = '/export/home/llecca/public/Gari/SLF123';

%% Probabilistic map generation
probMapGeneration(BASEDIR,'roisList','subjectList',...
    'testing','-MNISegment-linear','-linear-map');

%% Select for voxels with more than one ROI the ROI with higher number of 
%  subjects
%--------------------------------------------------------------------------
% 1)
% Binarize all the probabilistic maps so we can find where there are voxels
% that are labeled with more than one ROI

roi_thres = binarizeMaps(BASEDIR,'roisList','testing','-linear-map',...
    '-linear-binary','ROIs_sum');

% let's find voxels with are more than 1 ROI:
idx = find(roi_thres > 1);
[row, col, z] = ind2sub(size(roi_thres),idx);

% 2)
% Now let's obtain the sum of all the ROIs for all the sujects. This will 
% be useful to then count how many subjects are for each ROI in the 
% conflicting voxels
rois = countVoxels(BASEDIR,'roisList','subjectList',...
    '-MNISegment-linear-binary','countVoxels-testing');

% 3)
% Voting matrix generation
load(fullfile(BASEDIR,'countVoxels_ROIs_testing.mat'));
voting = votingMtrx(BASEDIR,'roisList',rois,row,col,z);

% now select for each voxel the ROI with maximum number of subjects: 
[~,idx] = max(voting,[],2);

% 4)
% Finally, clean the maps
cleanMaps(BASEDIR,'roisList','testing',idx,'-linear-map',...
    '-linear-thresholded-map',row,col,z);