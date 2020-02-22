function files = AFQ_mrtrix_build_files(fname_trunk,lmax,multishell)
% Builds a structure with the names of the files that the MRtrix commands
% will generate and need.
%
% files = mrtrix_build_files(fname_trunk,lmax)
%
% Franco Pestilli, Ariel Rokem, Bob Dougherty Stanford University
% GLU July.2016 added the T1 and tt5 file types
% GLU Jan.2019 added gmwmi and take it out with 5tt file

% Convert the raw dwi data to the mrtrix format: 
files.dwi = strcat(fname_trunk,'_dwi.mif');

% Extract the bo
files.b0 = strcat(fname_trunk,'_b0.mif');

% This file contains both bvecs and bvals, as per convention of mrtrix
files.b     = strcat(fname_trunk, '.b');

% If we have multishell data we will want the FA calculated with the shell
% closest to b1000
if multishell
    files.dwiSS = strcat(fname_trunk,'_dwiSS.mif');
    files.bSS   = strcat(fname_trunk, '.bSS');
end



% Convert the brain mask from mrDiffusion into a .mif file: 
files.brainmask         = strcat(fname_trunk,'_brainmask.mif'); 
files.brainmask_dilated = strcat(fname_trunk,'_brainmask_dilated.mif');
files.brainmask_eroded  = strcat(fname_trunk,'_brainmask_eroded.mif');

% Generate diffusion tensors:
files.dt = strcat(fname_trunk, '_dt.mif');

% Get the FA from the diffusion tensor estimates: 
files.fa = strcat(fname_trunk, '_fa.mif');

% Generate the eigenvectors, weighted by FA: 
files.ev = strcat(fname_trunk, '_ev.mif');

% Estimate the response function of single fibers: 
files.sf = strcat(fname_trunk, '_sf.mif');
files.response = strcat(fname_trunk, '_response.txt');

% Create a white-matter mask
files.wmMask         = strcat(fname_trunk, '_wmMask.mif');
files.wmMask_dilated = strcat(fname_trunk, '_wmMask_dilated.mif');

% Compute the CSD estimates: 
files.csd = strcat(fname_trunk, sprintf('_csd_lmax%i.mif',lmax)); 


% This was inside the if multishell before. Now create the 5tt files no matter what
    % Create tissue type segmentation to be used in multishell: 
        files.tt5   = strcat(fname_trunk, '_5tt.mif');
        files.gmwmi = strcat(fname_trunk, '_gmwmi.mif');
    
% The same with this. It is recommended now that the response function is calculated
% using dhollander function, so we need all
% if multishell
    % Create per tissue type response file
    files.wmResponse = strcat(fname_trunk, '_wmResponse.txt');
    files.gmResponse = strcat(fname_trunk, '_gmResponse.txt');
    files.csfResponse = strcat(fname_trunk, '_csfResponse.txt');
    % Compute the CSD estimates: 
    % (In this version we will use auto lmax, calculated by default by
    %                     dwi2fod, as per recommendation of dhollander)
    
    % files.wmCsd  = strcat(fname_trunk, sprintf('_wmCsd_lmax%i.mif',lmax)); 
    % files.gmCsd  = strcat(fname_trunk, sprintf('_gmCsd_lmax%i.mif',lmax)); 
    % files.csfCsd = strcat(fname_trunk, sprintf('_csfCsd_lmax%i.mif',lmax)); 
    
    files.wmCsd  = strcat(fname_trunk, '_wmCsd_autolmax.mif'); 
    files.gmCsd  = strcat(fname_trunk, '_gmCsd_autolmax.mif'); 
    files.csfCsd = strcat(fname_trunk, '_csfCsd_autolmax.mif'); 
    % RGB tissue signal contribution maps
    files.vf = strcat(fname_trunk, '_vf.mif');
    % dhollander voxel selection for QA
    files.voxels = strcat(fname_trunk, '_voxels.mif');
% end


end