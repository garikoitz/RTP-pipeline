function [status, results] = AFQ_mrtrix_maskfilter(files, ...
                                                   verbose, ...
                                                   mrtrixVersion)

%
% Calculate diffusion tensors.
%
% Parameters
% ----------
% in_file: The name of a dti file in .mif format
% out_file: The name of the resulting fa file in .mif format
% mask_file: The name of a mask file in .mif format (default to the entire
%             volume). 
% 
% Returns
% -------
% status: whether (0) or not (non-zero) the operation succeeded
% results: the results of the operation in stdout
%
% Notes
% -----
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
%
% Edited GLU 06.2016:
%        1.- Remove absolute paths
%        2.- Include mrTrix version


if notDefined('verbose')
    verbose = true;
end

if notDefined('bkgrnd')
    bkgrnd = false;
end





if mrtrixVersion == 2
    error('mrTrix version 2 is deprecated')
end
if mrtrixVersion == 3
    % Dilate brain mask
    % Adding option -npass we can further dilate it
    cmd_str1         = ['maskfilter -force ' files.brainmask  ' ' ...
                        'dilate - | mrthreshold -force -abs 0.5 - ' ...
                        files.brainmask_dilated];
    cmd_str2         = ['maskfilter -force ' files.brainmask  ' ' ...
                        '-npass 2 erode - | mrthreshold -force -abs 0.5 - ' ...
                        files.brainmask_eroded];

end


% Send it to mrtrix:
[status,results] = AFQ_mrtrix_cmd(cmd_str1, bkgrnd, verbose, mrtrixVersion);
[status,results] = AFQ_mrtrix_cmd(cmd_str2, bkgrnd, verbose, mrtrixVersion);

% in_file='/black/localhome/glerma/TESTDATA/AFQ/output/MareikeS03b/afq_20-Dec-2018_19h34m56s/dti96trilin/mrtrix/5tt_wm.mif';
% out_file='/black/localhome/glerma/TESTDATA/AFQ/output/MareikeS03b/afq_20-Dec-2018_19h34m56s/dti96trilin/mrtrix/5tt_wm_dilated.nii.gz';

