function [status, results] = AFQ_mrtrix_dwi2tensor(files, ...
                                                   multishell,...
                                                   bkgrnd, ...
                                                   verbose, ...
                                                   mrtrixVersion)

%
% Calculate diffusion tensors. 
%
% Parameters
% ----------
% in_file: The name of a diffusion file in .mif format
% out_file: The name of the resulting dti file in .mif format
% b_file: The name of a mrtrix format gradient vector file (default: 'encoding.b')
%
% Returns
% -------
% status: whether (0) or not (1) the operation succeeded 
% results: the results of the operation in the terminal
%
% Notes
% -----
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
% Latest edition: GLU 2019.02: if multishell use closest shell to 1000

if notDefined('verbose'); verbose = true; end
if notDefined('bkgrnd');  bkgrnd  = false;end
if mrtrixVersion ~= 3; error('Only mrTrix version 3 supported.');end


% This command generates  tensors: 
if multishell
    cmd_str = ['dwi2tensor -force ' ...
               '-grad ' files.bSS ' ' ...
               files.dwiSS ' ' ...
               files.dt];
else
    cmd_str = ['dwi2tensor -force ' ...
               '-grad ' files.b ' ' ...
               files.dwi ' ' ...
               files.dt];
end

% Send it to mrtrix: 
[status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose, mrtrixVersion); 