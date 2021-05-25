function rois = countVoxels(BASEDIR,roisList,subjectList,name2read,...
    countVoxels)

% This function summ all the subjects for each ROI and then locate the
% voxels with more than 1 ROI so we can see how many votes has each ROI
% -------------------------------------------------------------------------
% BASEDIR:      directory where is your data and the roisList.txt and
%               the subjectList.txt
% roisList:     .txt file name containing your ROIs names (do not include
%               the extension)
% subjectList:  .txt file containig your subjectList names (do not include
%               the extension)
% name2read:    it is the suffix of the name of the ROI to read, 
%               e.g. -MNISegment
% countVoxels:  name of the .mat file to save. It contains the struct
%               counting how many subjects are for each voxel and for each
%               ROI of the ROIsList
% -------------------------------------------------------------------------

fileID = fopen(fullfile(BASEDIR,strcat(roisList,'.txt')));
tlineRoi = fgetl(fileID);

roisum = 0;

while ischar(tlineRoi)
    fprintf('Loading ROI %s \n',tlineRoi)
    fileIDsub = fopen(fullfile(BASEDIR,strcat(subjectList,'.txt')));
    tlinesub = fgetl(fileIDsub);
    while ischar(tlinesub)
        fprintf('Subject %s \n',tlinesub)
        input_dir = fullfile(BASEDIR,tlinesub,'ROIs');
        roiFullPath = fullfile(input_dir,strcat(tlineRoi,name2read,...
            '.nii.gz'));
        % read roi for subject
        roi = niftiRead(roiFullPath);
        roisum = roi.data + roisum;
        tlinesub = fgetl(fileIDsub);
    end
    fclose(fileIDsub);
    rois.(tlineRoi)=roisum; roisum=0;
    tlineRoi = fgetl(fileID);
end
fclose(fileID);
save(fullfile(BASEDIR,strcat(countVoxels,'.mat')),'rois');

end