function roi = ANTS_CreateRoiFromMniNifti(roiNifti, warpNifti,refNifti, outfile)
% Transform a nifti roi saved in MNI space to an individuals brain

% Find the affine associated with the warp
f = strfind(warpNifti,'Inverse')-1;
affineXform = [warpNifti(1:f) 'Affine.txt'];

% If ROIs folder does not exist, create it
if ~exist(fileparts(outfile))
    mkdir(fileparts(outfile))
end

% Create the command to call ANTS
if exist('refNifti','var') && ~isempty(refNifti)
    % If a reference image was provided then reslice the roi to match this
    cmd = sprintf('WarpImageMultiTransform 3 %s %s -R %s --use-NN -i %s %s', roiNifti, outfile, refNifti, affineXform, warpNifti);
else
    cmd = sprintf('WarpImageMultiTransform 3 %s %s --use-NN -i %s %s', roiNifti, outfile, affineXform, warpNifti);
end
% excecute it
system(cmd);

% Convert nifti image to .mat roi
roi = dtiRoiFromNifti(outfile,[],prefix(prefix(outfile)),'mat');

% Clean the ROI
roi = dtiRoiClean(roi,3,{'fillHoles'});

% Save the ROI
dtiWriteRoi(roi,prefix(prefix(outfile)));

%% GLU: this is throwing me empty ROIs. Try with the code I was using in MINI
% This is an example call with the above code:
% WarpImageMultiTransform 3 
% /data/localhome/glerma/soft/vistasoft/mrDiffusion/templates/MNI_JHU_tracts_ROIs/ATR_roi1_L.nii.gz 
% /data/localhome/glerma/TESTDATA/AFQ/input/dtiInit_LTOZZI/dti150trilin/ROIs/ATR_roi1_L.nii.gz 
% --use-NN 
% -i /data/localhome/glerma/TESTDATA/AFQ/input/dtiInit_LTOZZI/dti150trilin/bin/b0Affine.txt 
% /data/localhome/glerma/TESTDATA/AFQ/input/dtiInit_LTOZZI/dti150trilin/bin/b0InverseWarp.nii.gz

% input = roiNifti;
% warp  = warpNifti;
%         cmd = str('antsApplyTransforms -d 3 '+
%                   '-i '+join(path2atlas,hemi+'-vols-1mm',nucleo)+' '+
%                   '-r '+fixedImage+' '+
%                   '-n Linear '+
%                   '-t '+join(opDir,'ants1Warp.nii.gz ')+
%                   '-t '+join(opDir,'ants0GenericAffine.mat ')+
%                   '-o '+join(opDir,hemi+'_'+nucleo))


% 

