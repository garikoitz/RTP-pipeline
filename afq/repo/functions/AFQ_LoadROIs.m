function [roi1, roi2] = AFQ_LoadROIs(fgNumber, subDir, ts)
% Load the two ROIs used to define a fiber group
% EDITED BY GLU: now it reads the ROIs from the tracts table, it is more flexible
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
% Update by Garikoitz Lerma March 2020
%    We cannot assume ROIs. They will need to exist as passed by the table. 
%    Deleting old code


% It doesn't matter which subject number because all we need is the roi
% name
RoiPara  = load(fullfile(subDir,'dt6.mat'));
fs_dir   = RoiPara.params.fs_dir;
roi_dir  = fullfile(fs_dir, 'ROIs');
roi1Name = char(fullfile(roi_dir, strcat(ts.roi1,ts.dilroi1,ts.extroi1)));
roi2Name = char(fullfile(roi_dir, strcat(ts.roi2,ts.dilroi2,ts.extroi2)));

if exist(roi1Name,'file') 
    roi1=dtiImportRoiFromNifti(roi1Name);
else
    error('ROI %s  does not exist for %s',roi1Name, subDir)   
end
if exist(roi2Name,'file')
    roi2=dtiImportRoiFromNifti(roi2Name);
else
    error('ROI %s  does not exist for %s',roi2Name, subDir)   
end

return
