function [status, results] = AFQ_mrtrix_brainmask(files, ...
                                                  bkgrnd, ...  
                                                  verbose, ...
                                                  mrtrixVersion)

%  GLU 02.2019


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
    % Create the b0 that we will copy as nifti to the /bin folder
    cmd_str = ['dwi2mask -force  ' ...
               '-grad ' files.b ' ' ...
                files.dwi ' ' ...
                files.brainmask];
   
end
% Send it to mrtrix:
[status,results] = AFQ_mrtrix_cmd(cmd_str,bkgrnd,verbose,mrtrixVersion);

