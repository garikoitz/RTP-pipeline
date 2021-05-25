function roi_thres = binarizeMaps(BASEDIR,roisList,folderoutput,...
    name2read,name2save,nameAllROis2save)

% This function binarize the probabilistic maps and save a nifti with all
% the ROIs added together, so we can find where there are voxels with
% more than 1 ROI, e.g., >1
% -------------------------------------------------------------------------
% BASEDIR:      directory where is your data and the roisList.txt and
%               the subjectList.txt
% roisLIst:     .txt file name containing your ROIs names (do not include
%               the extension)
% folderoutput: folder where are the probabilistic maps.
% name2read:    it is the suffix of the name of the ROI to read, 
%               e.g. -MNISegment
% name2save:    it is the suffix of the name of the map to save, e.g. -map
% nameAllROis2save: it is the suffix of the name of the map with all the 
%                   ROIs added together
% -------------------------------------------------------------------------
% OUTPUT:
% roi_thres: you can save this variable which contains the ROI with all the
% ROIs added together. this is useful to detect voxels contianing more than
% 1 ROI, e.g., in the data this voxels will have a value higher than 1
% -------------------------------------------------------------------------

fileID = fopen(fullfile(BASEDIR,strcat(roisList,'.txt')));
tlineRoi = fgetl(fileID);

roi_thres = 0;

% binarize the ROIs
while ischar(tlineRoi)
    fprintf('Loading ROI %s \n',tlineRoi)
    roiFullPath = fullfile(BASEDIR,folderoutput,strcat(tlineRoi,...
        name2read,'.nii.gz'));
    % read probabilistic map
    roi = niftiRead(roiFullPath);
    % set zero for values under zero
    roi.data(roi.data<0)=0;
    % set one for values above zero
    roi.data(roi.data>0)=1;
    roi_bin = roi.data;
    
    roi_thres = roi_thres + roi_bin;
    
    roi.cal_min = min(roi_bin(:));
    roi.cal_max = max(roi_bin(:));
    roi.fname = fullfile(BASEDIR,folderoutput,...
        strcat(tlineRoi,name2save,'.nii.gz'));
    writeFileNifti(roi)
    tlineRoi = fgetl(fileID);
end
roi.data = roi_thres;
roi.cal_min = min(roi_thres(:));
roi.cal_max = max(roi_thres(:));
roi.fname = fullfile(BASEDIR,folderoutput,...
    strcat(nameAllROis2save,'.nii.gz'));
writeFileNifti(roi);
fclose(fileID);

end

