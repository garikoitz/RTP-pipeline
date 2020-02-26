function [roi1, roi2] = AFQ_LoadROIs(fgNumber, subDir, ts)
% Load the two ROIs used to define a fiber group
% EDITED BY GLU: now it reads the ROIs from the tracts table, it is more flexible
% TODO: delete all the old code, first make it work with tables.
% 
% [roi1 roi2] = AFQ_LoadROIs(fgNumber,sub_dir. [afq])
%
% Inputs:
% fgNumber = The number of the fiber group. For example 1 is left ATR and
%            20 is right arcuate.
% sub_dir  = A path to the subjects directory
% afq      - An afq structure. This is optional but is needed if it is a
%            fiber group not in the mori groups
%
% Outputs:
% roi1     = First roi defining the tract
% roi2     = Second roi defining the tract
%
% Written by Jason D. Yeatman 1/31/2012

% If no afq structure was passed in it is fine for known ROIs
if false % ~exist('afq','var') || isempty(afq) && fgNumber <=20
    roi1Names={'ATR_roi1_L', 'ATR_roi1_R','CST_roi1_L','CST_roi1_R','CGC_roi1_L','CGC_roi1_R'...
        'HCC_roi1_L','HCC_roi1_R','FP_L','FA_L','IFO_roi1_L','IFO_roi1_R','ILF_roi1_L','ILF_roi1_R'...
        'SLF_roi1_L','SLF_roi1_R','UNC_roi1_L','UNC_roi1_R','SLF_roi1_L','SLF_roi1_R'};
    roi2Names={'ATR_roi2_L', 'ATR_roi2_R','CST_roi2_L','CST_roi2_R','CGC_roi2_L','CGC_roi2_R'...
        'HCC_roi2_L','HCC_roi2_R','FP_R','FA_R','IFO_roi2_L','IFO_roi2_R','ILF_roi2_L','ILF_roi2_R'...
        'SLF_roi2_L','SLF_roi2_R','UNC_roi2_L','UNC_roi2_R','SLFt_roi2_L','SLFt_roi2_R'};
    roiDir=fullfile(sub_dir,'ROIs');
    roi1path = fullfile(roiDir,[roi1Names{fgNumber} '.mat']);
    roi2path = fullfile(roiDir,[roi2Names{fgNumber} '.mat']);
    if exist(roi1path,'file') && exist(roi2path,'file')
        roi1=dtiReadRoi(roi1path);
        roi2=dtiReadRoi(roi2path);
    else
        error('ROI does not exist for %s',sub_dir)
    end
else
    % It doesn't matter which subject number because all we need is the roi
    % name
    RoiPara  = load(fullfile(subDir,'dt6.mat'));
    fs_dir   = RoiPara.params.fs_dir;
    roi_dir  = fullfile(fs_dir, 'ROIs');i
    roi1Name = char(fullfile(roi_dir, strcat(ts.roi1,ts.dilroi1,ts.extroi1)));
    roi2Name = char(fullfile(roi_dir, strcat(ts.roi2,ts.dilroi2,ts.extroi2)));





    if exist(roi1Name,'file') && exist(roi2Name,'file')
        % roi1=dtiReadRoi(roi1path);
        % roi2=dtiReadRoi(roi2path);
        roi1=dtiImportRoiFromNifti(roi1Name);
        roi2=dtiImportRoiFromNifti(roi2Name);

    else
        error('ROI does not exist for %s',sub_dir)
    end
end

return
