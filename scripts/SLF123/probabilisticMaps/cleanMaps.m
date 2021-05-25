function cleanMaps(BASEDIR,roisList,folderoutput,idx,name2read,...
    name2save,row,col,z)

% This function selects one ROI for the voxels with more than one ROI
% It selects the ROI with higher subjects for that voxel
% -------------------------------------------------------------------------
% BASEDIR:      directory where is your data and the roisList.txt and
%               the subjectList.txt
% roisList:     .txt file name containing your ROIs names (do not include
%               the extension)
% folderoutput: folder where to save the generated probabilistic maps. If
%               it does not exist, this function creates the folder.
% idx:          index refering the winner ROI for the voxel containing more
%               than one ROI
% name2read:    it is the suffix of the name of the ROI to read, 
%               e.g. -MNISegment
% name2save:    it is the suffix of the name of the map to save, e.g. -map
% row:          idx of the rows containing more than 1 voxel
% col:          idx of the columns contianing more than 1 voxel
% z:            idx of the heigh contaning more than 1 voxel
% -------------------------------------------------------------------------

fileID = fopen(fullfile(BASEDIR,strcat(roisList,'.txt')));
tlineRoi = fgetl(fileID);
count=0;

maps_dir = fullfile(BASEDIR,folderoutput);

while ischar(tlineRoi)
    fprintf('Loading ROI %s \n',tlineRoi);
    count=count+1;
    
    [a,~]=find(idx==count);
    a=setdiff(1:length(idx),a);
    
    roiFullPath = fullfile(BASEDIR,folderoutput,strcat(tlineRoi,...
        name2read,'.nii.gz'));
    roi = niftiRead(roiFullPath);
    for i=1:length(a)
        roi.data(row(a(i)),col(a(i)),z(a(i)))=0;
    end
    roi_bin=roi.data;
    roi.cal_min = min(roi_bin(:));
    roi.cal_max = max(roi_bin(:));
    roi.fname = fullfile(maps_dir,strcat(tlineRoi,name2save,'.nii.gz'));
    writeFileNifti(roi)
    
    tlineRoi = fgetl(fileID);
end

end

