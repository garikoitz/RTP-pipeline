function [status, results] = AFQ_mrtrix_5ttwm(tt5File, ...
                                              faFile, ...
                                              faMaskThresh, ...
                                              brainmask_eroded, ...
                                              wmMask, ...
                                              wmMask_dilated, ...
                                              verbose, ...
                                              mrtrixVersion, ...
                                              force)

% This is prepared to create the wmmask based on the anatomical 5tt file. 
% The problem is that if we do not have the reversed phase encoding file
% the alligment is not very good so the wmMask we are going to create it
% wont be perfect either...
% 
% Strategy;
% - Use FA > faMaskThresh to create a wmMask for seeding. 
% - Use eroded twice brainmask to clean the borders, it is usually messy
% - Use the wmMask obtained from 5tt file (anatomical) to create the mask
%   where the fibers cannot cross. 
% 
% 
% This is hte code to extract the wm only from the 5tt file
% 
% For now use the FA file to create the mask
%
%  GLU 01.2019


if notDefined('verbose')
    verbose = true;
end
if notDefined('bkgrnd')
    bkgrnd = false;
end
if notDefined('force'), force = false; end

if force
    f = '-force';
else
    f='';
end

if mrtrixVersion == 2
    error('mrTrix version 2 is deprecated')
end
if mrtrixVersion == 3
    % Create the fa based wmMask first, this will be for seeding
    % Using the eroded brain mask, clean the borders and save with the same name
    cmd_str1 = ['mrthreshold ' f ' -abs ' num2str(faMaskThresh) ' ' ...
                        faFile ' - | mrcalc ' f ' - ' ...
                        brainmask_eroded ' -mult ' wmMask];
        
    % Now create the anatomically based wmMask and dilate it, it will be used to 
    % contain the fibers, but as being it dilated it will encounter the GM
    
    % For now I am going to remove this, because I found that for some
    % datasets, the fsl segmentation doesn't work well...
    % And then the mask it is usually bigger than the wmMask, but it is not
    % perfect either, so I will use a bigger mask based in the previous
    % code now for testing and if it is ok it will stay for now. 
    % We will use ACT and the anatomicals when we have FS, otherwise the
    % analysis will be DWI based
        % cmd_str2 = ['mrconvert -coord 3 2 -axes 0,1,2 ' tt5File ' ' ...
        %          '- | maskfilter -force -npass 2 - dilate ' ...
        %           wmMask_dilated];
        
    cmd_str2 = ['mrthreshold ' f ' -abs ' num2str(faMaskThresh - .1) ' ' ...
                        faFile ' - | mrcalc ' f ' - ' ...
                        brainmask_eroded ' -mult ' wmMask_dilated];
        
    
end

% Send it to mrtrix:
[status,results] = AFQ_mrtrix_cmd(cmd_str1, bkgrnd, verbose,mrtrixVersion);
[status,results] = AFQ_mrtrix_cmd(cmd_str2, bkgrnd, verbose,mrtrixVersion);


