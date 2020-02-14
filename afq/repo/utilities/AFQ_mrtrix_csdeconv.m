function [status, results] = AFQ_mrtrix_csdeconv( files, ...
                                                  lmax, ...
                                                  bkgrnd, ...
                                                  verbose, ...
                                                  mrtrixVersion)
                                             
                                             
                                             
%function [status, results] = mrtrix_csdeconv(dwi_file, response_file, lmax, ...
%           out_file, grad_file, [mask_file=entire volume], [verbose=true])
%
% Fit the constrained spherical deconvolution model to dwi data 
%
% Parameters
% ----------
% dwi_file: The name of a dwi file in .mif format. 
% response_file: The name of a .txt fiber response function file. 
% lmax: The maximal harmonic order. 
% out_file: The resulting .mif file. 
% grad_file: a file containing gradient information in the mrtrix format. 
% mask_file: a .mif file containing a mask. Default: the entire volume. 
% verbose: Whether to display stdout to the command window. 
% 
% Returns
% -------
% status: whether (0) or not (1) the operation succeeded
% results: the results of the operation in stdout
%
% Notes
% -----
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
% 

if notDefined('verbose')
    verbose = true; 
end

if notDefined('bkgrnd')
    bkgrnd = false;
end

mrtrixVersion = 3

% According to D Tournier, this is what he recommends for single shell
% If you’re using single-shell CSD [...] a better option would be to switch to 2-tissue CSD:
%        estimate responses for GM, WM, CSF using dwi2response dhollander
%        use dwi2fod msmt_csd using the WM & CSF responses only.
% This is what I am implementing here.
% In another thread dhollander recommends not including the lmax either,
% because the algo takes care of it. I removed the lmax info from the filename

cmd_str = ['dwi2fod msmt_csd ' ...
                     files.dwi ' ' ...
                     files.wmResponse  ' '  files.wmCsd ' '  ...
                     files.csfResponse ' '  files.csfCsd ' '  ...
                     '-mask ' files.wmMask_dilated ' ' ...
                     '-grad ' files.b ' ' ...
                     '-force'];                              
[status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);



% Esto es lo que he lanzado en command line
% dwi2fod msmt_csd -mask data_aligned_trilin_noMEC_brainmask.mif data_aligned_trilin_noMEC_dwi.mif data_aligned_trilin_noMEC_wm.txt data_aligned_trilin_noMEC_wm.mif data_aligned_trilin_noMEC_gm.txt data_aligned_trilin_noMEC_gm.mif data_aligned_trilin_noMEC_csf.txt data_aligned_trilin_noMEC_csf.mif -grad data_aligned_trilin_noMEC.b




