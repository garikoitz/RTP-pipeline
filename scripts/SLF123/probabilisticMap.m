% load vistasoft functions
addpath(genpath('/export/home/llecca/llecca/Scripts/SLF123/vistasoft'));

% define subjects database
BASEDIR = '/export/home/llecca/llecca/Scripts/SLF123';

% read the subjects to convert
fileID = fopen(strcat(BASEDIR,'/roisList.txt'));
tlineRoi = fgetl(fileID);

roiprob = 0;
counter = 0;

% create probabilisticMaps folder if does not exist
maps_dir = fullfile(BASEDIR,'probabilisticMaps');
if ~exist(maps_dir,'dir')
    mkdir(maps_dir);
end

while ischar(tlineRoi)
    fprintf('Loading ROI %s for all subjects list \n',tlineRoi)
    
    fileIDsub = fopen(strcat(BASEDIR,'/subjectList.txt'));
    tlinesub = fgetl(fileIDsub);
    
    while ischar(tlinesub)
        input_dir = strcat(BASEDIR,'/',tlinesub);
        roiFullPath = fullfile(input_dir,'ROIs',strcat(tlineRoi,...
            '-MNISegment.nii.gz'));
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
    roi.fname = fullfile(maps_dir,strcat(tlineRoi,'-map.nii.gz'));
    writeFileNifti(roi)
    roiprob = 0;
    counter = 0;
    tlineRoi = fgetl(fileID);
end
fclose(fileID);
