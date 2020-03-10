function files = AFQ_mrtrixInit(dt6, ...
                                lmax, ...
                                mrtrix_folder, ...
                                mrtrixVersion, ...
                                multishell, ...
                                tool, ...
                                faMaskThresh, ...
                                force)
% function files = AFQ_mrtrixInit(dt6, lmax, mrtrix_folder)
% 
% Initialize an mrtrix CSD analysis
%
% This fucntion computes all the files needed to use mrtrix_track. 
%
% Parameters
% ----------
% dt6: string, full-path to an mrInit-generated dt6 file.
% T1nii: path to the acpc-ed T1w nii used at the beginning. 
% lmax: The maximal harmonic order to fit in the spherical deconvolution (d
%       model. Must be an even integer. This input determines the
%       flexibility  of the resulting model fit (higher values correspond
%       to more flexible models), but also determines the number of
%       parameters that need to be fit. The number of dw directions
%       acquired should be larger than the number of parameters required.
% mrtrix_folder; Name of the output folder
%
% Notes
% -----
% This performs the following operations:
%
% 1. Convert the raw dwi file into .mif format
% 2. Convert the bvecs, bvals into .b format
% 3. Convert the brain-mask to .mif format 
% 4. Fit DTI and calculate FA and EV images
% 5. Estimate the response function for single fibers, based on voxels with
%    FA > 0.7
% 6. Fit the CSD model. 
% 7. Convert the white-matter mask to .mif format. 
% 
% For details: 
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
% 
% Update GLU 2018.06:
% When trying to dockerize it, it is not working as some of the paths of
% the dt6 come hardcoded. Edit it to make everything relative paths 


if notDefined('mrtrix_folder'), mrtrix_folder = 'mrtrix'; end
if notDefined('lmax'), lmax = 4; end
if notDefined('tool'), tool = 'freesurfer'; end
if notDefined('multishell'), multishell = false; end
if notDefined('mrtrixVersion'), mrtrixVersion = 3; end
if notDefined('faMaskThresh'), faMaskThresh = 0.3; end
if notDefined('force'), force = true; end

% Loading the dt file containing all the paths to the fiels we need.
dt_info = load(dt6);
disp('Just load dt6, and the dt6 = dt_info.files is:')
dt_info.files
%

% Check if the dt6 coming from the dtiinit GEAR and the afq GEAR are
% running in the same space. Otherwise fix the paths
dt6Parts = split(dt_info.files.alignedDwRaw, filesep);
mrtrixFolderParts  = split(mrtrix_folder, filesep);
% If the first element is different, assume that they come from different
% spaces, recreate the folders and filenames
if ~strcmp(dt6Parts{2}, mrtrixFolderParts{2})
    warning('DtiInit and AFQ have been processed in different environments')
end


% Obtain the session name. This is usually the zip name if it has not
% been edited. 
SessionDir = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-2)), filesep)
% This is where the dt6 is located
AnalysisDir = fileparts(dt6);
AnalysisDir2 = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-1)), filesep);
% Not very useful sanity check:
if ~strcmp(AnalysisDir, AnalysisDir2); error('Check all the path mess for mrtrix files...'); end

% Strip the file names out of the dt6 strings. 
% dwRawFile = dt_info.files.alignedDwRaw;
% dwRawFile = fullfile(dt_info.params.rawDataDir, strcat(dt_info.params.rawDataFile,'.gz'));
% dwRawFile = fullfile(SessionDir, strcat(dt_info.params.rawDataFile,'.gz'));
% The first option is the good one, but if I did this change I want to maintain it. Nevertheless, it cannot come from dt_info.params
[p,f,e] = fileparts(dt_info.files.alignedDwRaw);
dwRawFile = fullfile(SessionDir,[f e]); 


% This line removes the extension of the file (.nii.gz) and mainaints de path
fname_trunk = dwRawFile(1:strfind(dwRawFile,'.')-1);
% With this code we can separate the rest
[pathDwRawFile, fnameDwRawFile] = fileparts(fname_trunk);

% In the mrtrix_folder argument we already have the path to the mrtrix
% folder
mrtrix_dir = mrtrix_folder;

% Assuming in 'session' we want the subject_name/dmri64 or whatever
session = pathDwRawFile;
% session = dt_info.params.rawDataDir;

% And in fname_trunk we want the whole path and the beginning of the
% filename
fname_trunk = [mrtrix_folder filesep fnameDwRawFile]; 

if ~exist(mrtrix_dir, 'dir')
    mkdir(mrtrix_dir)
end

% Build the mrtrix file names.
files = AFQ_mrtrix_build_files(fname_trunk, lmax, multishell);

% Check wich processes were already computed and which ones need to be done.
computed = mrtrix_check_processes(files);

% Convert the raw dwi data to the mrtrix format: 
if (~computed.('dwi')) || (~computed.('b')) 
    % dwRawFile
    % files.dwi
	[bvecPath bvecName bvecExt] = fileparts(dt_info.files.alignedDwBvecs);
    [bvalPath bvalName bvalExt] = fileparts(dt_info.files.alignedDwBvals);
    bvecs = fullfile(session, [bvecName bvecExt]);
    bvals = fullfile(session, [bvalName bvalExt]);
    cmd   = sprintf('mrconvert -force -fslgrad %s %s -stride 1,2,3,4  -export_grad_mrtrix %s %s %s', ...
					 bvecs, bvals, files.b, dwRawFile, files.dwi);
    [status,results] = AFQ_mrtrix_cmd(cmd, 0,0,3);
end


% Create the b0: do it always, so pass multishell variable as false
if (~computed.('b0'))
    AFQ_mrtrix_extract(files,false,0,0,mrtrixVersion);
end
% Only when multishell extract the shell closest to 1000, we will use this
% files only to calculate the dt and fa
if multishell
    if (~computed.('dwiSS') || (~computed.('bSS')))
        AFQ_mrtrix_extract(files,true,0,0,mrtrixVersion);
    end
end

% Convert the brain mask from mrDiffusion into a .mif file.
% We are not creating hte brainmask anymore, create it here
if (~computed.('brainmask'))
    AFQ_mrtrix_brainmask(files,0,0,mrtrixVersion);
end

% Dilate and erode the brainmask
if (~computed.('brainmask_dilated'))  || (~computed.('brainmask_eroded'))
  brainMaskFile        = fullfile(session, dt_info.files.brainMask); 
  AFQ_mrtrix_maskfilter(files, ...
                        false, ...
                        mrtrixVersion);
end

% Generate diffusion tensors:
if (~computed.('dt'))
  AFQ_mrtrix_dwi2tensor(files,...
                        multishell, ...
                        0, ...
                        0, ...
                        mrtrixVersion);
end

% Get the FA from the diffusion tensor estimates: 
if (~computed.('fa'))
  AFQ_mrtrix_tensor2FA(files, ...
                       0, ...
                       0, ...
                       mrtrixVersion);
end

% Generate the eigenvectors, weighted by FA: 
if  (~computed.('ev'))
  AFQ_mrtrix_tensor2vector(files.dt, files.ev, files.fa,0,mrtrixVersion);
end

% Estimate the response function of single fibers: 
if (~computed.('wmResponse')) & (~computed.('gmResponse')) & (~computed.('csfResponse'))
  AFQ_mrtrix_response(files, ... 
                      false, ... % this is show_figure
                      [], ...    % bckground
                      lmax, ...  % lmax: not necessary, they recommend to leave it blank, but for future proof leave it here
                      false, ... % verbose
                      mrtrixVersion, ...
                      force) 
end

% Create the 5tt file from the T1 / aparc+aseg
if (~computed.('tt5')) || (~computed.('gmwmi'))
    inputFile = [];
    switch tool
        case {'fsl'}        
            inputFile = fullfile(session, dt_info.files.t1);
            if ~(exist(inputFile, 'file') == 2)
                error(['Cannot find T1, please copy it to ' session]);
            end    
            % Find and aseg file and better if it is an aparc+aseg one. 
            % Select the first aparc if there are several *aseg* files.
            % It can take mgz or nii or mif
        case {'freesurfer'}
            % Add it to the dt_info for future versions
            asegFiles = dir(fullfile(session,'*aseg*'));
            for ii = 1:length(asegFiles)
               if length(strfind(asegFiles(ii).name, 'aseg')) > 0
                   inputFile = fullfile(session, asegFiles(ii).name);
               end
               if length(strfind(asegFiles(ii).name, 'aparc')) > 0
                   inputFile = fullfile(session, asegFiles(ii).name);
               end
            end
            if ~(exist(inputFile, 'file') == 2)
                disp(['inputFile = ' inputFile]); 
                warning(['Cannot find aseg file, please copy it to ' session ...
                         'if there is a t1 file fsl will be used']);
                tool = 'fsl';
                inputFile = fullfile(session, dt_info.files.t1);
                if ~(exist(inputFile, 'file') == 2)
                    error(['Cannot find T1, please copy it to ' session]);
                end 
            end
        otherwise
            error(sprintf('The tool %s has not been implemented, use freesurfer or fsl.', tool))
    end
    fprintf('Running the 5ttgen code with the %s tool and %s\n', tool, inputFile)
    AFQ_mrtrix_5ttgen(inputFile, ...
                      files.tt5, ...
                      files.gmwmi, ...
                      0, ...
                      0, ...
                      mrtrixVersion,...
                      tool, ...
                      force);
end

% Create a white-matter mask
if (~computed.('wmMask')) || (~computed.('wmMask_dilated'))
    AFQ_mrtrix_5ttwm(files.tt5, ...
                      files.fa, ...
                      faMaskThresh, ...
                      files.brainmask_eroded, ...
                      files.wmMask, ...
                      files.wmMask_dilated, ...
                      false, ...
                      mrtrixVersion, ...
                      force)

end

% Calculate the FOD using CSD
% This is separated in two different calls now, but with the recent changes
% in the mrtrix code, they should be in the same one...
if multishell
    % Create per tissue response function estimation
    % Not using the other response function, we already have the masks
    % 2019: updating and moving it to the response function to have it all
    % together, following new mrTrix recommendations
     
    % Compute the CSD estimates: 
    if ~computed.('wmCsd')
      disp('The following step takes a while (a few hours)');                                  
      AFQ_mrtrix_csdeconv_msmt(files, ...
                               lmax, ...
                               0, ...
                               0, ...
                               mrtrixVersion)
    end
    
    % RGB tissue signal contribution maps
    if ~computed.('vf')
         % mrconvert -coord 3 0 wm.mif - | mrcat csf.mif gm.mif - vf.mif
        cmd_str = ['mrconvert -coord 3 0 ' files.wmCsd ' - | ' ...
                   'mrcat ' files.gmCsd ' ' files.csfCsd ' - ' files.vf];           
        AFQ_mrtrix_cmd(cmd_str, 0, 0,mrtrixVersion);
    end
    
else
    % Compute the CSD estimates: 
    if (~computed.('wmCsd'))  
      disp('The following step takes a while (a few hours)'); 
      disp('Using the default CSD with dilated brainmask, recommended for ACT, it only takes a little bit longer'); 
      AFQ_mrtrix_csdeconv(files, ...
                          lmax, ...
                          false,... % Verbose
                          false, ...
                          mrtrixVersion)
    end
end

end  % End function



