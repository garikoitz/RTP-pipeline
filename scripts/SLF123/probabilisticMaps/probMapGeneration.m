function probMapGeneration(BASEDIR,roisList,subjectList,folderoutput,...
    name2read,name2save)

% This function generates probabilistc maps from all the subjects defined 
% in the subjectList.txt. It will be a map for each ROI denfined in the
% ROIsList.txt
% -------------------------------------------------------------------------
% BASEDIR:      directory where is your data and the roisList.txt and
%               the subjectList.txt
% roisLIst:     .txt file name containing your ROIs names (do not include
%               the extension)
% subjectList:  .txt file containig your subjectList names (do not include
%               the extension)
% folderoutput: folder where to save the generated probabilistic maps. If
%               it does not exist, this function creates the folder.
% name2read:    it is the suffix of the name of the ROI to read, 
%               e.g. -MNISegment
% name2save:    it is the suffix of the name of the map to save, e.g. -map
% -------------------------------------------------------------------------

fileID = fopen(fullfile(BASEDIR,strcat(roisList,'.txt')));
tlineRoi = fgetl(fileID);

roiprob = 0;
counter = 0;

% create probabilisticMaps folder if does not exist
maps_dir = fullfile(BASEDIR,folderoutput);
if ~exist(maps_dir,'dir')
    mkdir(maps_dir);
end

while ischar(tlineRoi)
    fprintf('Loading ROI %s for all subjects list \n',tlineRoi)
    
    fileIDsub = fopen(fullfile(BASEDIR,strcat(subjectList,'.txt')));
    tlinesub = fgetl(fileIDsub);
    
    while ischar(tlinesub)
        input_dir = fullfile(BASEDIR,tlinesub,'ROIs');
        roiFullPath = fullfile(input_dir,strcat(tlineRoi,...
            name2read,'.nii.gz'));
        % read roi for subject
        roi = niftiRead(roiFullPath);
        % set zero for values under zero
        roi.data(roi.data<0)=0;
        roi.data = roi.data ./ max(roi.data(:));
        % add subject data to the roi
        roiprob = roi.data+roiprob;
        % subject counter
        counter = counter+1; disp(counter);
        tlinesub = fgetl(fileIDsub);
    end
    fclose(fileIDsub); 
    roiprob = (roiprob ./ counter).*100;
    % change img data
    roi.data = round(roiprob);
    roi.cal_min = min(roiprob(:));
    roi.cal_max = max(roiprob(:));
    roi.fname = fullfile(maps_dir,strcat(tlineRoi,name2save,'.nii.gz'));
    writeFileNifti(roi)
    roiprob = 0;
    counter = 0;
    tlineRoi = fgetl(fileID);
end
fclose(fileID);

end

