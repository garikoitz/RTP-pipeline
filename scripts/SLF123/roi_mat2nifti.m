% ROI.mat to ROI.nii.gz conversion

% load vistasoft functions
addpath(genpath('/export/home/llecca/llecca/Scripts/SLF123/vistasoft'));

% define subjects database
BASEDIR = '/bcbl/home/public/Gari/SLF123';

% read the subjects to convert
fileID = fopen(strcat(BASEDIR,'/subSesList.txt'));
% below we start a loop
tline = fgetl(fileID);
while ischar(tline)
    fprintf('Working on: %s \n',tline)
    input_dir=strcat(BASEDIR,'/',tline);
    checkfile = fullfile(input_dir,'/dt6/bin/','b0.nii.gz');
    if exist(checkfile,'file')
        fprintf('%s exists \n',checkfile)
        img = niftiRead(checkfile);
        rois = dir(fullfile(input_dir, 'ROIs', '*.mat'));
        for df=1:length(rois)
            roiFullPath = fullfile(input_dir, 'ROIs',rois(df).name);
            roi         = dtiReadRoi(roiFullPath);
            coords      = roi.coords;
            % Convert vertex acpc coords to img coords
            imgCoords  = mrAnatXformCoords(img.qto_ijk, coords);
            % Get coords for the unique voxels
            imgCoords = unique(ceil(imgCoords),'rows');
            % Make a 3D image
            roiData = zeros(img.dim);
            roiData(sub2ind(img.dim, imgCoords(:,1), imgCoords(:,2), imgCoords(:,3))) = 1;
            % Change img data
            img.data = roiData;
            img.cal_min = min(roiData(:));
            img.cal_max = max(roiData(:));
            % Write the nifti file
            [~,roiNameWoExt] = fileparts(roiFullPath);
            img.fname = fullfile(fileparts(roiFullPath), [roiNameWoExt,'.nii.gz']); 
            writeFileNifti(img); 
        end
    else
        fprintf('%s does not exist \n',checkfile)
    end
    tline = fgetl(fileID);
end
fclose(fileID);

