function [status, results] = AFQ_mrtrix_5ttgen(input_file, ...
                                                  tt5_filename, ...
                                                  gmwmi_filename, ...
                                                  bkgrnd, ...
                                                  verbose, ...
                                                  mrtrixVersion,...
                                                  tool, ...
                                                  force) 
%
% Convert a nifti image to a mrtrix .mif image.
% Do this: 5ttgen fsl T1.mif 5tt.mif
% To do: explore more options of the software
% Example
% AFQ_mrtrix_5ttgen('T1.mif', '5tt.mif') 
% 
% Sthg strange is going on with cluster tmp folders
% I create a /tmp folder in $HOME for every subject. 
% The script deletes all folders, but it wil maintain the $HOME/tmp 

if notDefined('bkgrnd'),  bkgrnd  = false;end
if notDefined('verbose'), verbose = true; end
if notDefined('mrtrixVersion'), mrtrixVersion = 3; end
if notDefined('tool'), tool = 'fsl'; end
if notDefined('force'), force = false; end

if force
    f = '-force';
else
    f='';
end

% There were problems in the cluster with the /tmp directory, and there
% have been reports with the same problems. Just in case use a tmp dir in
% the home folder. Everything will be deleted automatically after use. 
% If you want to maintain the tmp folders for visualization add -nocleanup
baseDir = fileparts(tt5_filename);
tmpDir = fullfile(baseDir,'tmp');
if ~exist(tmpDir, 'dir'), mkdir(tmpDir), end



% fsl or freesurfer can be selected
if strcmp(tool, 'fsl')
    cmd_str = ['5ttgen fsl ' f '  ' ...
               input_file ' ' ...
               tt5_filename ' ' ...
               '-nocrop -scratch ' tmpDir];
else
    lutPath = flutlocation();
    cmd_str = ['5ttgen freesurfer ' f ' ' ...
               '-lut ' lutPath ' ' ...
               input_file   ' ' ...
               tt5_filename ' ' ...
               '-nocrop -scratch ' tmpDir];
    
end

[status,results] = AFQ_mrtrix_cmd(cmd_str, ...
                                  bkgrnd, ...
                                  verbose, ...
                                  mrtrixVersion); 
                              

% Create the seed with this 5tt2gmwmi

cmd_str   = ['5tt2gmwmi ' f ' ' tt5_filename  ' ' gmwmi_filename ];
AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose, mrtrixVersion)


end

