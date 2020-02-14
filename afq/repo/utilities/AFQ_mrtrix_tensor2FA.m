function [status, results] = AFQ_mrtrix_tensor2FA(files, ...
                                                  bkgrnd, ...
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
% Edited GLU 02.2019:
%        1.- Use eroded to avoid FA in skull


if notDefined('verbose'); verbose = true; end
if notDefined('bkgrnd');  bkgrnd  = false;end
if mrtrixVersion ~= 3; error('Only mrTrix version 3 supported.');end



cmd_str = ['tensor2metric -force ' ...
           '-fa - ' files.dt ' | ' ...
               'mrcalc - ' files.brainmask_eroded ' ' ...
               '-mult ' files.fa];



% Send it to mrtrix:
[status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);