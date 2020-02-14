function [status,results] = AFQ_mrtrix_response(files, ...
                                                 show_figure, ...
                                                 bkgrnd,  ...
                                                 lmax, ...
                                                 verbose, ...
                                                 mrtrixVersion)
% Calculate the fiber response function utilized by MRtrix for the spherical
% deconvolution analysis.
%
% INPUTS
%     files     - The structure with all the filenames created in
%                 mrtrix_build_files
%   show_figure - Optional. Whether to show a figure of the response function
%                 profile (default: true)
%       verbose - Outputs the opretions being performted in the MatLab
%                 prompt.
% 
% OUTPUTS
%        status - whether (0) or not (1) the operation succeeded
%       results - the results of the operation in the terminal
%
% NOTES
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
%
% Franco Pestilli, Bob Dougherty and Ariel Rokem, Stanford University 2013
% Garikoitz Lerma-Usabiaga, BCBL 2016
% Garikoitz Lerma-Usabiaga, Stanford University 2019


if notDefined('verbose'),             verbose = true;end
if notDefined('bkgrnd'),               bkgrnd = false;end
if notDefined('lmax'),                   lmax = 6; end
if notDefined('mrtrixVersion'), mrtrixVersion = 3;end


% D TOurnier recommends here using the MS even for single shell
% http://community.mrtrix.org/t/spurious-streamlines-difference-between-tckgen-versions/1705/4
% estimate responses for GM, WM, CSF using dwi2response dhollander
% use dwi2fod msmt_csd using the WM & CSF responses only.
% Here, it is recommended not specify lmax or mask
% http://community.mrtrix.org/t/dwi2response-dhollander-questions/1095

% It seems that regardless of the anatomical registration, dhollander
% performs better than the msmt 5tt option, so, in conclusion, per every
% dataset, bein MS or not, having rpe or not, we are going to use
% dhollander to calculate the response function

% Build the generic cmd_str that will always be applied
cmd_str = ['dwi2response dhollander -force ' ...
                          files.dwi ' ' ...
                          files.wmResponse  ' ' ...
                          files.gmResponse  ' ' ...
                          files.csfResponse ' ' ...
                          '-grad ' files.b  ' ' ...
                          '-voxels ' files.voxels];

[status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose, mrtrixVersion);    

end


% This is the function I used to run directly in AFQ_mrTrixInit.m function. 
% if (~computed.('wmResponse')) && (mrtrixVersion > 2)
%     cmd_str = ['dwi2response msmt_5tt ' ...
%                 files.dwi ' ' files.tt5 ' ' ...
%                 files.wmResponse ' ' files.gmResponse ' ' files.csfResponse ' ' ...
%                 '-grad ' files.b];
% 
%     AFQ_mrtrix_cmd(cmd_str, 0, 0,mrtrixVersion);
% end



% Example called used for MS data in 2016, to be updated
% dwi2response msmt_5tt data_aligned_trilin_noMEC_dwi.mif data_aligned_trilin_noMEC_5tt.mif data_aligned_trilin_noMEC_wm.txt data_aligned_trilin_noMEC_gm.txt data_aligned_trilin_noMEC_csf.txt -grad data_aligned_trilin_noMEC.b



