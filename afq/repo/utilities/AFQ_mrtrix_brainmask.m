function [status, results] = AFQ_mrtrix_brainmask(files, ...
                                                  bkgrnd, ...  
                                                  verbose, ...
                                                  mrtrixVersion, ...
                                                  fs_brainmask)

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
    %{
    cmd_str = ['dwi2mask -force  ' ...
               '-grad ' files.b ' ' ...
                files.dwi ' ' ...
                files.brainmask];
     %}   
   % dwi2mask is problematic, so we are going to use the one from fs
   % the rest of the steps will be the same but with a differnt file. 
   cmd = ['mrconvert -force -stride 1,2,3,4  ' fs_brainmask ' - | ' ...
          'mrtransform -force -template ' files.b0 ' - - | ' ...
          'mrthreshold -force -abs 0.5 - ' files.brainmask];
       
end
% Send it to mrtrix:
[status,results] = AFQ_mrtrix_cmd(cmd,bkgrnd,verbose,mrtrixVersion);
% Check
if ~isfile(files.brainmask);error('%s could not be created',files.brainmask);end
