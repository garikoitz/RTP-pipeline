function voting = votingMtrx(BASEDIR,roisList,rois,row,col,z)

% This function generates the voting array (i.e., the array counting how
% many subjects for each ROI are for the voxels containing more than 1 ROI)
% -------------------------------------------------------------------------
% BASEDIR:  directory where is your data and the roisList.txt and the 
%           subjectList.txt
% roisLIst: .txt file name containing your ROIs names (do not include
%           the extension)
% rois:     struct containing all the ROIs added together to see the
%           voxels with more than 1 ROI
% row:      idx of the rows containing more than 1 voxel
% col:      idx of the columns contianing more than 1 voxel
% z:        idx of the heigh contaning more than 1 voxel
% -------------------------------------------------------------------------

fileID = fopen(fullfile(BASEDIR,strcat(roisList,'.txt')));
tlineRoi = fgetl(fileID);

voting = zeros(length(z),size(struct2table(rois),2));
count=0;

while ischar(tlineRoi)
    fprintf('Loading ROI %s \n',tlineRoi);
    roi = rois.(tlineRoi);
    count=count+1;
    for t=1:length(z)
        voting(t,count) = roi(row(t),col(t),z(t));
    end
    tlineRoi = fgetl(fileID);
end
fclose(fileID);

end

