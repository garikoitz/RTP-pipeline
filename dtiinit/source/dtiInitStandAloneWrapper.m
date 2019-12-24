function dtiInitStandAloneWrapper(json)
% 
% dtiInitStandAloneWrapper(json)
%
% Read a JSON object, a JSON file, or a directory containing a json file
% and run dtiInit inside of a docker container (vistalab/dtiinit). 
% 
% 
% INPUTS:
%       json - a JSON string, a JSON file, or a directory containing a json
%              file, in the following format (Note that 'input_dir' and
%              'output_dir' are the only REQUIRED inputs)
%
% OUTPUTS: 
%       A docker run produces a zip file containing all of the outputs
%       from the algorithm. The name of the output zip file is:
%           'dtiInit[date-time].zip'
%
%
% JSON SCHEMA:
%       Below is an example JSON file with the defaults show for 'params'.
%       See dtiInitParams.m for more info regarding params. Note that
%       "input_dir" and "output_dir" are required and must be in the
%       context of the container. 
%{
            { 
                "input_dir": "/input",
                "output_dir": "/output",
                "dwi_file": "",
                "bvec_file": "",
                "bval_file": "",
                "t1_file": "",
                "aparcaseg_file": "",
                "params":  
                    {
                        "bvalue": "",
                        "gradDirsCode": "",
                        "clobber": 0,
                        "dt6BaseName": "",
                        "flipLrApFlag": 0,
                        "numBootStrapSamples": 500,
                        "fitMethod": "ls",
                        "nStep": 50,
                        "eddyCorrect": 0,
                        "excludeVols": "",
                        "bsplineInterpFlag": 0,
                        "phaseEncodeDir": "",
                        "dwOutMm": [2, 2, 2],
                        "rotateBvecsWithRx": 0,
                        "rotateBvecsWithCanXform": 0,
                        "bvecsFile": "",
                        "bvalsFile": "",
                        "noiseCalcMethod": "b0",
                        "outDir": "/output/"
                    }
            }
%}
% 
% REQUIRED INPUTS:
%       'input_dir' and 'output_dir' are the only required inputs.  
% 
% 
% HELP: 
%       If 'help', '-h', '--help', or nothing (nargin==0), is passed in
%       this help will be displayed.
% 
% 
% USAGE:
%       Pass in a JSON file, a JSON text string, or a path to a directory
%       containing a JSON file to the docker container to initiate a
%       dtiInit processing run (see INPUT section for JSON schema):
% 
%       % Using a JSON file
%        docker run --rm -ti -v `pwd`/input:/input -v `pwd`/output:/output vistalab/dtiinit /input/dtiInit.json
% 
%       % Using a JSON string
%        docker run --rm -ti -v `pwd`/input:/input -v `pwd`/output:/output vistalab/dtiinit '{"input_dir":"/input", "output_dir": "/output"}'
% 
%       % Using a directory (in the container), containing a JSON (.json)
%        docker run --rm -ti -v `pwd`/input:/input -v `pwd`/output:/output vistalab/dtiinit /input/
% 
% Use compile.sh for compiling
% Use this command to launch in matlab
%{

jsonPath = '/data/localhome/glerma/soft/RTP-pipeline/dtiinit/source';
% jsonPath = '~/soft/rtp-pipeline/dtiinit/source';

  % dtiInitStandAloneWrapper(fullfile(jsonPath,'dtiInit.json'));
  % dtiInitStandAloneWrapper(fullfile(jsonPath,'dtiInit_with_aparcAseg_4ltozziMS.json'));
  % dtiInitStandAloneWrapper(fullfile(jsonPath,'dtiInit_with_aparcAseg_4HCPMS.json'));
  % dtiInitStandAloneWrapper(fullfile(jsonPath,'defining.json'));
  %   dtiInitStandAloneWrapper(fullfile(jsonPath,'HCPdep.json'));
      dtiInitStandAloneWrapper(fullfile(jsonPath,'HCPdepMBP.json'));
%}
% Use this command to run the docker in the directory
% 
% (C) Vista Lab, Stanford University, 2015
% 
% TODO:
%   - Reproducibility
%   - Remote download of input files
%   - Json Structure validation
%   - Zip the results.

%% Initial checks

% If nothing was passed in, display help and return
if nargin == 0
    help_file = '/opt/help.txt';
    if exist(help_file, 'file')
        system(['cat ', help_file]);
    else
        help(mfilename);
    end
    return
end

% Assume the user wanted to see the help, and show it
if ischar(json) 
    if strcmpi(json, 'help') || strcmpi(json, '-help') || strcmpi(json, '-h') || strcmpi(json, '--help')
        help(mfilename);
    end
end


%% Parse the JSON file or object

if exist(json, 'file') == 2
    J = loadjson(json);
elseif exist(json, 'dir') == 7
    jsonFile = dir(fullfile(json, '*.json'));
    jsonFile = fullfile(json, jsonFile.name);
    disp(jsonFile);
    if ~isempty(jsonFile)
        J = loadjson(jsonFile);
    else
        error('No JSON file could be found');
    end
elseif ~isempty(json) && ischar(json)
    try
        J = loadjson(json);
    catch ME
        disp(ME.message); 
        return
    end
else
    error('Could not find/parse the json file/structure');
end


%% Check the json ojbect for required fields

required = {'input_dir', 'output_dir'};
err = false;

for r = 1:numel(required)
    if ~isfield(J, required{r})
        err = true;
        fprintf('%s not found in JSON object!\n', required{r});
    elseif ~exist(J.(required{r}), 'dir')
        fprintf('%s Does not exist!\n', required{r});
        err = true;
    end
end

% If there was a problem, return
if err 
    error('Exiting! There was a problem with the inputs. Please check input_dir and output_dir!');
end

% Create an output subfolder for the outputs 
outputSubFolder = ['dtiInit_', strrep(strrep(datestr(now),' ', '_'),':','-')];
J.output_dir = fullfile(J.output_dir, outputSubFolder);
mkdir(J.output_dir);


%% Get a list of diffusion files from the input directory

dw = getDwiFilesStruct(J.input_dir);
dw = dw{1}; % For this case, limit to the first set found

if ~isfield(J, 'dwi_file') || ~exist(J.dwi_file, 'file')
    J.dwi_file = dw.nifti;
end
if ~isfield(J, 'bvec_file') || ~exist(J.bvec_file, 'file')
    J.bvec_file = dw.bvec;
end
if ~isfield(J, 'bval_file') || ~exist(J.bval_file, 'file')
    J.bval_file = dw.bval;
end


%% T1 File (Paths for templates are container specific)

if ~isfield(J, 't1_file') || ~exist(J.t1_file, 'file')
    template_t1 = '/templates/MNI_EPI.nii.gz'; 
    J.t1_file = template_t1;
end


%% Initialize diffusion parameters

dwParams            = dtiInitParams;
dwParams.outDir     = J.output_dir;
dwParams.bvecsFile  = J.bvec_file;
dwParams.bvalsFile  = J.bval_file;
dwParams.bvalue     = dw.bvalue;

% disp('This is dwParams: ', dwParams)

%% Update the diffusion params from the JSON object

if isfield(J, 'params')
    param_names = fieldnames(J.params);
    for f = 1:numel(param_names)
        if isfield(dwParams,param_names{f}) && ~isempty(J.params.(param_names{f}))
            dwParams.(param_names{f}) = J.params.(param_names{f});
        end
    end
else
    disp('Using default dtiInit params')
end

disp(J)
disp(J.params)
disp(dwParams)


%% Validate the JSON structure against the JSON schema
% GLU REMOVED
% fprintf('Validating JSON input against schema... ');
% dtiInitStandAloneValidateJson(J);
% fprintf('Success.\n');


% Validate that the bval values are normalized
% If not, write again the file normalized so that it can be read downstream
% Determine shell
bvals = dlmread(J.bval_file);
roundedBval  = 100 * round(bvals/100);
paramsShells = unique(roundedBval);
if 0 == min(paramsShells)
    paramsShells = paramsShells(paramsShells ~= 0);
    numShells    = length(paramsShells);
else
    error('It seems that this file have no b0. Check it please.')
end

% Write the files back
warning('The bVals were normalized.')
dlmwrite(J.bval_file, roundedBval, 'delimiter',' ');


%% Run dtiInit
% From 3.0.5 onwards I am forking dtiInit and giving it less
% functionalities. I will remove the tensor fitting and do it with mrTrix,
% it is much faster and the rest of the code relies in it anyways. 
% First iteration, stop dtiInit doing it and later on mrtrixInit will paste
% the fa to the bin folder. 
% In 3.1.1 it was doing nothing already
% In 3.1.2 I am removing the call altogether to make this thing simpler
%        AFQ_dtiInit(J.dwi_file, J.t1_file, dwParams);
% Plans for 3.2.0: remove afq-browser and dtiInit altogether


% Here dtiInit was called, assign variables
dwRawFileName = J.dwi_file;
t1FileName    = J.t1_file;

% I. Load the diffusion data, set up parameters and directories structure
% Load the difusion data
disp('Loading preprocessed (use rtp-preproc or already preprocessed) data...');
dwRaw = niftiRead(dwRawFileName);

% By default all processed nifti's will be at the same resolution as the
% dwi data
% if notDefined('dwParams'); 
%   dwParams         = dtiInitParams; 
% we are going to ignore whatever we pass in the dtinit params. 
% remove them from manifest too in next version
  dwParams.dwOutMm = dwRaw.pixdim(1:3);
% end 

% Initialize the structure containing all directory info and file names
dwDir      = dtiInitDir(dwRawFileName,dwParams);
outBaseDir = dwDir.outBaseDir;
fprintf('Dims = [%d %d %d %d] \nData Dir = %s \n', size(dwRaw.data), dwDir.dataDir);
fprintf('Output Dir = %s \n', dwDir.subjectDir);


% II. Select the anatomy file

% Check for the case that the user wants to align to MNI instead of T1.
if exist('t1FileName','var') && strcmpi(t1FileName,'MNI')
    t1FileName = fullfile(mrDiffusionDir,'templates','MNI_EPI.nii.gz');
    disp('The MNI EPI template will be used for alignment.');
end

if notDefined('t1FileName') || ~exist(t1FileName,'file')
    t1FileName = mrvSelectFile('r',{'*.nii.gz';'*.*'},'Select T1 nifti file');
    if isempty(t1FileName); disp('dtiInit canceled by user.'); return; end
end
fprintf('t1FileName = %s;\n', t1FileName);


% New in 3.1.2: check the file is RAS If not, convert to RAS Be careful, bvecs
% needs to be changed as well accordingly Instead of changing bvecs manually, we
% will let mrtrix take care of it Convert files to mif, add the bvecs and bvals
% as part of the file, do the conversion and then output the fsl type bvecs
% again. This we can maintain coherence in the whole process. Therefore, if
% there is a - in the strides or the order is not 123, we need to convert (with
% freesurfer it is easier because it gives you RAS or PIR or LAS or whatever but
% it is not installed in the Docker container)

% Check it in these files
checkfiles = {'dwi_file','t1_file','aparcaseg_file'};
for nc=1:length(checkfiles)
    fname = J.(checkfiles{nc});
    [c2r,orientation] = rtp_convert2RAScheck(fname);
    if ~c2r
        fprintf('%s has orientation %s (RAS), no transformation required.\n',fname,orientation)
    else
        fprintf('%s has orientation %s, converting to 123(4) (RAS) \n',fname,orientation)
        [p,f,e] = fileparts(fname);
        % Check that it is .nii.gz
        if strcmp(e,'.gz')
            [~,f] = fileparts(f);
            e = '.nii.gz';
        elseif strcmp(e,'.nii')
            % do nothing
        else
            error('File extension should be .nii.gz or .nii')
        end
        % If it is diffusion, we need to take care of gradients
        switch checkfiles{nc}
            case 'dwi_file'
                %{
                % 1: convert it to mif with integraded gradient info
                mrtname    = fullfile(p,[f '.mif']);
                cmd = sprintf('mrconvert -quiet -force -fslgrad %s %s %s %s', J.bvec_file, J.bval_file, fname, mrtname);
                s = system(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not convert to mif');end
                % 2: convert it to RAS strides
                mrtRASname = fullfile(p,[f 'ras.mif']);
                cmd = sprintf('mrconvert -quiet -force -stride 1,2,3,4 %s %s', mrtname, mrtRASname);
                s = system(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not convert to RAS');end
                % 3: convert back to nifti and export fsl gradients to go back
                %    to original situation
                fslRASname = fullfile(p,[f 'ras.nii.gz']);
                fslRASbvec = fullfile(p,[f 'ras.bvec']);
                fslRASbval = fullfile(p,[f 'ras.bval']);
                cmd = sprintf('mrconvert -quiet -force -export_grad_fsl %s %s %s %s', fslRASbvec, fslRASbval, mrtRASname, fslRASname);
                s = system(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not convert back to nifti RAS');end
                % 4: remove the .mif files
                % cmd = sprintf('rm -f %s',mrtname);
                % s = system(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not remove mif file');end
                % cmd = sprintf('rm -f %s',mrtRASname);
                % s = system(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not remove mif file');end
                % 5: change filenames to continue with processing normally
                J.bvec_file        = fslRASbvec;
                J.bval_file        = fslRASbval;
                J.(checkfiles{nc}) = fslRASname;
                dwParams.bvecsFile = J.bvec_file;
                dwParams.bvalsFile = J.bval_file;
                %}
                
                
                
                
                
                
                % 1: Do the strides to RAS conversion
                fslRASname = fullfile(p,[f 'ras.nii.gz']);
                fslRASbvec = fullfile(p,[f 'ras.bvec']);
                fslRASbval = fullfile(p,[f 'ras.bval']);
                cmd = sprintf('mrconvert -quiet -force -fslgrad %s %s -stride 1,2,3,4  -export_grad_fsl %s %s %s %s', ...
                              J.bvec_file, J.bval_file,fslRASbvec,fslRASbval,fname, fslRASname);
                
                s = AFQ_mrtrix_cmd(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not change the strides to RAS');end
                % 2: change filenames to continue with processing normally
                J.bvec_file        = fslRASbvec;
                J.bval_file        = fslRASbval;
                J.(checkfiles{nc}) = fslRASname;
                dwRawFileName      = J.(checkfiles{nc});
                dwParams.bvecsFile = J.bvec_file;
                dwParams.bvalsFile = J.bval_file;
            otherwise
                % 1: Do the strides to RAS conversion
                fRASname = fullfile(p,[f 'ras.nii.gz']);
                cmd = sprintf('mrconvert -quiet -force -stride 1,2,3 %s %s', fname, fRASname);
                s = AFQ_mrtrix_cmd(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not change the strides to RAS');end
                % 2: change filenames to continue with processing normally
                J.(checkfiles{nc}) = fRASname;
        end
    end
end












% XV. Name the folder that will contain the dt6.mat file
% GLU REVISE: 
% If the user passed in a full path to dt6BaseName and outDir ... if
% they're different the dt6.mat file will be saved to dt6BaseName while the
% other data will be saved to outDir. See dtiInitDir for the fix.
% Remove the nUniqueDirs thing as well, why should it be in the folder name
if isempty(dwParams.dt6BaseName) 
    % nUniqueDirs from dtiBootGetPermMatrix
    % dwParams.dt6BaseName = fullfile(dwDir.subjectDir,sprintf('dti%02d',nUniqueDirs));
    dwParams.dt6BaseName = fullfile(dwDir.subjectDir,'dticsd');
    %if ~dwParams.bsplineInterpFlag 
        % Using trilinear interpolation 
     %   dwParams.dt6BaseName = [dwParams.dt6BaseName 'trilin'];
    %end
else
    if isempty(fileparts(dwParams.dt6BaseName)) 
        dwParams.dt6BaseName = fullfile(dwDir.subjectDir,dwParams.dt6BaseName);
    end
end


% XVII. Build the dt6.files field and append it to dt6.mat
% GLU: This was the old version before I commented XVI
% Need to handle the case where there is more than one dt6 file. 
% for dd = 1:numel(dt6FileName)
%     dtiInitDt6Files(dt6FileName{dd},dwDir,t1FileName);
% end

% GLU: the new one
outBaseName = dwParams.dt6BaseName;
dt6FileName = fullfile(outBaseName, 'dt6.mat');
binDirName  = fullfile(outBaseName, 'bin');
if(~exist(outBaseName,'dir'));mkdir(outBaseName);end
if(~exist(binDirName,'dir')) ;mkdir(binDirName);end
if(~exist('adcUnits','var')); adcUnits = ''; end

params.buildDate = datestr(now,'yyyy-mm-dd HH:MM');
l = license('inuse');
params.buildId = sprintf('%s on Matlab R%s (%s)',l(1).user,version('-release'),computer);
if(ischar(dwRawFileName)); [dataDir,rawDataFileName] = fileparts(dwRawFileName);  % dwRawAligned);
else                      [dataDir,rawDataFileName] = fileparts(dwRawFileName.fname);end  % dwRawAligned.fname); end
endparams.rawDataDir = dataDir;
params.rawDataFile   = rawDataFileName;
% We assume that the raw data file is a directory inside the 'subject' directory.
params.subDir = fileparts(dataDir);



% Some day I will try to understand why they are doing this...
[fullParentDir, binDir] = fileparts(binDirName);
[ppBinDir, pBinDir] = fileparts(fullParentDir);
pBinDir = fullfile(pBinDir,binDir);


% Now decide which ones of this files I will create with mrTrix and which
% ones I will leave uncreated
files.b0        = fullfile(pBinDir,'b0.nii.gz');
files.brainMask = fullfile(pBinDir,'brainMask.nii.gz');
files.wmMask    = fullfile(pBinDir,'wmMask.nii.gz');
% files.wmProb    = fullfile(pBinDir,'wmProb.nii.gz');
files.tensors   = fullfile(pBinDir,'tensors.nii.gz');
files.fa        = fullfile(pBinDir,'fa.nii.gz');
% files.vecRgb    = fullfile(pBinDir,'vectorRGB.nii.gz');
% files.faStd     = fullfile(pBinDir,'faStd.nii.gz');
% files.mdStd     = fullfile(pBinDir,'mdStd.nii.gz');
% files.pddDisp   = fullfile(pBinDir,'pddDispersion.nii.gz');

% This is new, we are going to copy the input data as output data and call it alligned in dt6
copyfile(dwRawFileName, dwParams.outDir);
[~,fname,ext] = fileparts(dwRawFileName);
files.alignedDwRaw   = fullfile(dwParams.outDir, [fname ext]);
dwDir.dwAlignedRawFile = files.alignedDwRaw;

copyfile(dwParams.bvecsFile, dwParams.outDir);
[~,fname,ext] = fileparts(dwParams.bvecsFile);
files.alignedDwBvecs = fullfile(dwParams.outDir,[fname ext]); 
dwDir.alignedBvecsFile = files.alignedDwBvecs;

copyfile(dwParams.bvalsFile, dwParams.outDir);
[~,fname,ext] = fileparts(dwParams.bvalsFile);
files.alignedDwBvals = fullfile(dwParams.outDir,[fname ext]);
dwDir.alignedBvalsFile = files.alignedDwBvals;


save(dt6FileName,'adcUnits','params','files');
dtiInitDt6Files(dt6FileName,dwDir,t1FileName);


% XX. Save out parameters, svn revision info, etc. for future reference

dtiInitLog(dwParams,dwDir);



% Exit operations
disp('***************** IMPORTANT *****************')
disp(sprintf('Copying the following files from dtiInit to AFQ: %s and %s',J.t1_file,J.aparcaseg_file))
disp('***************** IMPORTANT *****************')
copyfile(J.t1_file, J.output_dir)
copyfile(J.aparcaseg_file, J.output_dir)


%% Permissions

fileattrib(J.output_dir,'+w +x', 'o'); 


%% Compress the outputs

fprintf('Compressing output [%s]... ', J.output_dir);
cd(mrvDirup(J.output_dir));
zip([outputSubFolder, '.zip'], J.output_dir);
fprintf('Done.\n');


%% Remove uncompressed output files

rmdir(J.output_dir, 's');


%% TODO: REPRODUCIBILITY!


return 


