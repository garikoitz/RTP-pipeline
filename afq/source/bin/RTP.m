function RTP(jsonargs)
%
% AFQ_StandAlone_QMR(jsonargs)
%
% INPUT ARGUMENTS: Note: must be a json string.
%
%       input_dir:  Directory containing dt6.mat
%       out_name:   Name for the resulting afq file
%       output_dir: Location to save the results
%       params:     Key-Value pair of params for AFQ_Create
%       metadata:   Key_Value pairs for analysis.
%                   Defaults are:
%                           age = '';
%                           sex = '';
%                           ndirs = 30;
%                           bvalue = 1000;
%                           age_comp = false;
%{
% EXAMPLE USAGE:
%       jsonargs = '{"input_dir": "/home/lmperry/AFQ_docker/mrDiffusion_sampleData/dti40", "out_name": "afq.mat", "output_dir": "/black/lmperry/AFQ_docker/mrDiffusion_sampleData/dti40/AFQ3" }'
%       jsonargs = '{"input_dir" : "/data/localhome/glerma/TESTDATA/AFQ/input/dtiInit222/dti90trilin", "output_dir": "/data/localhome/glerma/TESTDATA/AFQ/output/withDtiinit222_mrtrix","params"    :"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json"}'%       jsonargs = '{"input_dir" : "/data/localhome/glerma/TESTDATA/AFQ/input/dtiInit111/dti90trilin", "output_dir": "/data/localhome/glerma/TESTDATA/AFQ/output/withDtiinit111_mrtrix","params"    :"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json"}'
%       jsonargs = '{"input_dir" : "/data/localhome/glerma/TESTDATA/AFQ/input/dtiInit222/dti90trilin", "output_dir": "/data/localhome/glerma/TESTDATA/AFQ/output/withDtiinit222_mrtrix","params"    :"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json"}'
%       jsonargs = '{"input_dir" : "/data/localhome/glerma/TESTDATA/AFQ/input/MareikeS13/dti96trilin", "output_dir": "/data/localhome/glerma/TESTDATA/AFQ/output/MareikeS13","params"    :"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json"}'
%       jsonargs = '{"input_dir" : "/data/localhome/glerma/TESTDATA/AFQ/input/MareikeS13act/dti96trilin", "output_dir": "/data/localhome/glerma/TESTDATA/AFQ/output/MareikeS13act","params"    :"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json"}'
jsonargs = ['{"input_dir" :' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/MareikeS13/dti96trilin",' ...
            '"output_dir": ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/output/MareikeS13", ' ...
            '"params"    : ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json"}']
jsonargs = ['{"input_dir" :' ...
            '"/Volumes/users/glerma/TESTDATA/AFQ/input/MareikeS13/dti96trilin",' ...
            '"output_dir": ' ...
            '"/Volumes/users/glerma/TESTDATA/AFQ/output/MareikeS13", ' ...
            '"params"    : ' ...
            '"/Volumes/users/glerma/TESTDATA/AFQ/input/config_parsed.json"}']

jsonargs = ['{"input_dir" :' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/dtiInit_LTOZZI/dti150trilin",' ...
            '"output_dir": ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/output/ltozzi", ' ...
            '"params"    : ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json"}']
% Testing the defining white matter tracts from Kurt Schilling
jsonargs = ['{"input_dir" :' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/dtiInit_defining/dti90trilin/",' ...
            '"output_dir": ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/output/defining", ' ...
            '"params"    : ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed_defining.json"}']     
% Testing the defining white matter tracts from Bertsolari
jsonargs = ['{"input_dir" :' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/dtiInit_bertso/dti64trilin/",' ...
            '"output_dir": ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/output/bertso", ' ...
            '"params"    : ' ...
            '"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed_bertso.json"}'] 
AFQ_StandAlone_QMR(jsonargs);

jsonargs = ['{"input_dir" :' ...
            '"/Users/glerma/soft/rtp-pipeline/local/DTIINIT/output/dtiInit_24-Dec-2019_07-55-23/dticsd/",' ...
            '"output_dir": ' ...
            '"/Users/glerma/soft/rtp-pipeline/local/AFQ/output", ' ...
            '"params"    : ' ...
            '"/Users/glerma/soft/rtp-pipeline/local/AFQ/input/config_parsed.json"}'] 



jsonargs = ['{"params"    : ' ...
            '"/black/localhome/glerma/soft/RTP-pipeline/example_output.json"}'] 
RTP(jsonargs);
%}
%
%#ok<*AGROW>


%% Begin

disp('Starting AFQ...');

% Initialize args
input_dir  = [];
out_name   = [];
output_dir = [];
params     = [];
% metadata   = [];

%% Handle jsonargs
disp('This is the json string to be read by loadjson:')
disp(jsonargs)



if exist('jsonargs', 'var') && ~isempty(jsonargs);
    P = loadjson(jsonargs);

end

%% Parse the params and setup the AFQ structure
%{
if ~isempty(params)
    if ischar(params)
        P = loadjson(params);
    else
        P = params;
    end
end
%}
P
%% Configure inputs and defaults
input_dir = P.input_dir;
output_dir = P.output_dir;
if notDefined('input_dir')
    if exist('/input', 'dir')
        input_dir = '/input';
    else
        error('An input directory was not specified.');
    end
end
sub_dirs{1} = input_dir;

if notDefined('output_dir')
    if exist('/output', 'dir')
        output_dir = '/output';
    else
        error('An output directory was not specified.');
    end
end
output_dir = fullfile(output_dir, 'AFQ');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
% update ROI directory
P.params.fs_dir = P.fs_dir;

% Just one group here
sub_group = ones(numel(sub_dirs),1);





%% See if I can read the templates
disp('_____ CHECK IF IT CAN READ TEMPLATES  _______')
% Get the AFQ base directory
AFQbase = AFQ_directories
% Template directory
tdir = fullfile(AFQbase,'templates','labelMaps')
% Path to the template
Tpath = fullfile(tdir,'MNI_AAL_AndMore.nii.gz')
% Load the template
if exist(Tpath,'file')
    Timg = readFileNifti(Tpath)
else
    error('Cannot read %s', Tpath)
end



%% CHECK INPUTS AND OPTIONALS
% LMX: until now I think it should have worked. Here we need to be sure first
%      that we can read the input files correctly. 
% The path and the variable names below might seems arbitrtraty, but we need to
% try to maintain them the same names because we don't know where they might
% fail down the line. There are old functions that still rely on the dt6.mat
% thing to know where the files are. If we decide to change it (what we snould
% do maybe in a later version) we should do it very throughly. 

% The next section is new, the RAS conversion, but you'll see that I am
% maintaining the ugly paths and old system. There is no J.files structure at
% this point, so we can create it here. 



% I AM GOING TO COPY HERE THE RELEVANT CODE COMING FROM DTI INIT

% File and folder checks
% LMX: not going to copy here, bad code. Check the files exist and the folder exist
% otherwise create them. You can check the code, but it basically assummes that
% the filenames are the same otherwise it breaks, etc. I prefer to be explicit.
% Give the actual correct names or break, I don't want to work but with the
% wrong files and not knowing. You can check dtiInitStandAloneWrapper line 180,
% but I would just ignore it.


% Templates: 
% This is the default template is using right now, it is very bad... add it for
% now but we will substitute it
% Variable names might need adjusting

%% copy input files to output/AFQ/
if(~exist(fullfile(P.output_dir,'AFQ'),'dir'));mkdir(fullfile(P.output_dir,'AFQ'));end
copyfile(fullfile(P.anat_dir,'t1.nii.gz'), fullfile(P.output_dir,'AFQ/'));
copyfile(fullfile(P.bvec_dir,'dwi.bvecs'), fullfile(P.output_dir,'AFQ/'));
copyfile(fullfile(P.bval_dir,'dwi.bvals'), fullfile(P.output_dir,'AFQ/'));
copyfile(fullfile(P.nifti_dir,'dwi.nii.gz'), fullfile(P.output_dir,'AFQ/'));
copyfile(fullfile(P.fs_dir, 'aparc+aseg.nii.gz'), fullfile(P.output_dir,'AFQ/'));
J.aparcaseg_file = fullfile(P.output_dir,'AFQ','aparc+aseg.nii.gz');
J.t1_file = fullfile(P.output_dir,'AFQ', 't1.nii.gz');
J.bvec_file = fullfile(P.output_dir,'AFQ', 'dwi.bvecs');
J.bval_file = fullfile(P.output_dir,'AFQ', 'dwi.bvals');
J.dwi_file = fullfile(P.output_dir,'AFQ', 'dwi.nii.gz');

J.input_dir = fullfile(P.output_dir,'AFQ');
J.output_dir = P.output_dir;
sub_dirs{1} = J.input_dir;
if ~isfield(J, 't1_file') || ~exist(J.t1_file, 'file')
    template_t1 = '/templates/MNI_EPI.nii.gz'; 
    J.t1_file = template_t1;
end

%% Initialize diffusion parameters
% LMX: maintain this structure for now, we'll see if we need it later
dwParams            = dtiInitParams;
dwParams.outDir     = J.output_dir;
dwParams.bvecsFile  = J.bvec_file;
dwParams.bvalsFile  = J.bval_file;
%dwParams.bvalue     = dw.bvalue;

%% Update the diffusion params from the JSON object
% LMX: ths was reading the dtiinit params, that we are not reading now but that
% they are in manifest.json. I think we don't need them, maintain this commented
% here and we will decide what to delete and what to add to the afq json params.
% 
% if isfield(J, 'params')
%     param_names = fieldnames(J.params);
%     for f = 1:numel(param_names)
%         if isfield(dwParams,param_names{f}) && ~isempty(J.params.(param_names{f}))
%             dwParams.(param_names{f}) = J.params.(param_names{f});
%         end
%     end
% else
%     disp('Using default dtiInit params')
% end
% % Print what is being passed to do sanity check when running in Docker mode
% disp(J)
% disp(J.params)
% disp(dwParams)

%% Validate that the bval values are normalized
% LMX II added this because downstream it was required to be normaliized. Maybe
% we can make it optional later. 
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
% We are in 4.0.0 now, that we already removed dtiInit and afq-browser.
%   As you see, from 3.1.2 I copied all the functionalities into this file, so
%   everything should be done here. Be careful with filenames and paths and the
%   like, before this was a different container, now both are ni the same
%   container, and some of the stuff I was doing was to be sure the right files
%   were passed from container to container. 


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
% LMX: WE SHOULD NEVER DO THIS! but, once I had a subjet without T1 so I used
% it so...
if exist('t1FileName','var') && strcmpi(t1FileName,'MNI')
    t1FileName = fullfile(mrDiffusionDir,'templates','MNI_EPI.nii.gz');
    disp('The MNI EPI template will be used for alignment.');
end

if notDefined('t1FileName') || ~exist(t1FileName,'file')
    t1FileName = mrvSelectFile('r',{'*.nii.gz';'*.*'},'Select T1 nifti file');
    if isempty(t1FileName); disp('dtiInit canceled by user.'); return; end
end
fprintf('t1FileName = %s;\n', t1FileName);


% (lmx: see that I am moving the old code with notes and everything, once we
% make it work, we'll clearn everything, but just in case it helps us debugging
% I am going to maintain it for now)

% XV. Name the folder that will contain the dt6.mat file
% If the user passed in a full path to dt6BaseName and outDir ... if
% they're different the dt6.mat file will be saved to dt6BaseName while the
% other data will be saved to outDir. See dtiInitDir for the fix.
% I removed the nUniqueDirs thing as well, why should it be in the folder name
if isempty(dwParams.dt6BaseName) 
    % nUniqueDirs from dtiBootGetPermMatrix
    % dwParams.dt6BaseName = fullfile(dwDir.subjectDir,sprintf('dti%02d',nUniqueDirs));
    dwParams.dt6BaseName = fullfile(dwDir.dataDir,'dticsd');
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

J.params = P.params;
J.params.input_dir = J.input_dir;
J.params.output_dir = J.output_dir;
params = J.params;
params.buildDate = datestr(now,'yyyy-mm-dd HH:MM');
l = license('inuse');
params.buildId = sprintf('%s on Matlab R%s (%s)',l(1).user,version('-release'),computer);
if(ischar(dwRawFileName)); [dataDir,rawDataFileName] = fileparts(dwRawFileName);  % dwRawAligned);
else                      [dataDir,rawDataFileName] = fileparts(dwRawFileName.fname);end  % dwRawAligned.fname); end
endparams.rawDataDir = dataDir;
params.rawDataFile   = rawDataFileName;
% We assume that the raw data file is a directory inside the 'subject' directory.
%params.subDir = fileparts(dataDir);
params.subDir = dataDir;


% Some day I will try to understand why they are doing this...
[fullParentDir, binDir] = fileparts(binDirName);
[ppBinDir, pBinDir] = fileparts(fullParentDir);
pBinDir = fullfile(pBinDir,binDir);


% Now decide which ones of this files I will create with mrTrix and which
% ones I will leave uncreated

% LMX: MAINTAIN THIS HERE. Here I just created the path names to the files that
% were created in dtiinit. So I just created the paths in dtiinit and then ni
% afq using mrtrix I created the files themselves. Rigth now we should change th
%e order of things, but it shuoold continue working

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

% LMX: Check this one, if we have been leaving the files in the correct place,
% this should not be necessary, check. 

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

% LMX: this is required. we are saving the dt6.mat file with all the required
% variables and file names, required for the rest of the thing
save(dt6FileName,'adcUnits','params','files');
dtiInitDt6Files(dt6FileName,dwDir,t1FileName);
copyfile(dt6FileName, J.input_dir);

% XX. Save out parameters, svn revision info, etc. for future reference

dtiInitLog(dwParams,dwDir);



% Exit operations (lmx: check what is necessary here, probably nothing, as we don't need to pass files between conotainers now)
disp('***************** IMPORTANT *****************')
disp(sprintf('Copying the following files from dtiInit to AFQ: %s and %s',J.t1_file,J.aparcaseg_file))
disp('***************** IMPORTANT *****************')
copyfile(J.t1_file, J.output_dir);
copyfile(J.aparcaseg_file, J.output_dir);

% (HERE  IT ENDS WHAT IT WAS IN DTI INIT)

%% CHECK IF INPUT IS RAS  (THIS WAS IN AFQ CONTAINER)
% New in 3.1.2: check the file is RAS If not, convert to RAS Be careful, bvecs
% needs to be changed as well accordingly Instead of changing bvecs manually, we
% will let mrtrix take care of it Convert files to mif, add the bvecs and bvals
% as part of the file, do the conversion and then output the fsl type bvecs
% again. This we can maintain coherence in the whole process. Therefore, if
% there is a - in the strides or the order is not 123, we need to convert (with
% freesurfer it is easier because it gives you RAS or PIR or LAS or whatever but
% it is not installed in the Docker container)
% This was originially done in the dtiinit wrapper but mrtrix did not work
% there... It was anoother step in removing the whoole dtiinit thing


% (lmx: now you can adapt this, you already now what are the files and where they are)
% First thing we do once we now where we are and we have the dt6.mat file, is
% checking that the files are RAS. If they are coming from rtp-preproc they will
% be, but it is not always the case. So we need to check. 

% Read the input file names and convert them
input_dir = J.input_dir;
basedir = split(input_dir,filesep);
basedir = strcat(basedir(1:end-1));
basedir = fullfile('/',basedir{:});

fprintf('This is input_dir: %s\n', input_dir)
fprintf('This is basedir: %s\n',   basedir)

J       = load(fullfile(input_dir,'dt6.mat'));
disp('This are the contents of dt6.mat')
J

% DWI file
[p,f,e]= fileparts(J.files.alignedDwRaw);
if ~strcmp(p,basedir); J.files.alignedDwRaw = fullfile(basedir,[f e]); end

% BVEC file
[p,f,e]= fileparts(J.files.alignedDwBvecs);
if ~strcmp(p,basedir); J.files.alignedDwBvecs = fullfile(basedir,[f e]); end

% BVAL file
[p,f,e]= fileparts(J.files.alignedDwBvals);
if ~strcmp(p,basedir); J.files.alignedDwBvals = fullfile(basedir,[f e]); end


J.files.t1path    = fullfile(basedir,J.files.t1);
fprintf('This is the absolute path to the t1: %s\n', J.files.t1path)
if exist(J.files.t1path,'file')
    fprintf('T1 file %s Exists. \n', J.files.t1path)
else
    error('Cannot find %s', J.files.t1path)
end
J.files.aparcaseg = fullfile(basedir,J.files.t1);
% Solve the aparc+aseg case
asegFiles = dir(fullfile(basedir,'*aseg*'));
for ii = 1:length(asegFiles)
    if length(strfind(asegFiles(ii).name, 'aseg')) > 0
        J.files.aparcaseg = fullfile(basedir, asegFiles(ii).name);
    end
    if length(strfind(asegFiles(ii).name, 'aparc')) > 0
        J.files.aparcaseg = fullfile(basedir, asegFiles(ii).name);
    end
end
if ~(exist(J.files.aparcaseg, 'file') == 2)
    disp(['inputFile = ' J.files.aparcaseg]);
    warning(['Cannot find aseg file, please copy it to ' basedir]);
    
    J.files.aparcaseg = J.files.t1path;
    if ~(exist(J.files.aparcaseg, 'file') == 2)
        error(['Cannot find T1, please copy it to ' basedir]);
    end
end

% Check it in these files
checkfiles = {'alignedDwRaw','t1path','aparcaseg'};
for nc=1:length(checkfiles)
    fname = J.files.(checkfiles{nc});
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
            case 'alignedDwRaw'
                % 1: Do the strides to RAS conversion
                fslRASname = fullfile(p,[f 'ras.nii.gz']);
                fslRASbvec = fullfile(p,[f 'ras.bvec']);
                fslRASbval = fullfile(p,[f 'ras.bval']);
                cmd = sprintf('mrconvert -quiet -force -fslgrad %s %s -stride 1,2,3,4  -export_grad_fsl %s %s %s %s', ...
                              J.files.alignedDwBvecs, J.files.alignedDwBvals,fslRASbvec,fslRASbval,fname, fslRASname);
                
                s = AFQ_mrtrix_cmd(cmd); if s~=0;error('[RTP] Could not change the strides to RAS');end
                % 2: change filenames to continue with processing normally
                J.files.alignedDwBvecs = fslRASbvec;
                J.files.alignedDwBvals = fslRASbval;
                J.files.alignedDwRaw   = fslRASname;
                % dwRawFileName      = J.(checkfiles{nc});
                % dwParams.bvecsFile = J.bvec_file;
                % dwParams.bvalsFile = J.bval_file;
            otherwise
                % 1: Do the strides to RAS conversion
                fRASname = fullfile(p,[f 'ras.nii.gz']);
                cmd = sprintf('mrconvert -quiet -force -stride 1,2,3 %s %s', ...
                               fname, fRASname);
                s = AFQ_mrtrix_cmd(cmd); if s~=0;error('[dtiInitStandAloneWrapper] Could not change the strides to RAS');end
                % 2: change filenames to continue with processing normally
                J.files.(checkfiles{nc}) = fRASname;
                
        end
    end
end
% Change t1 back
[p,f,e] = fileparts(J.files.t1path);
J.files.t1 = [f,e];
% Save the new dt6 file 
adcUnits = '';
params   = J.params;
files    = J.files;
save(fullfile(input_dir,'dt6.mat'),'adcUnits','params','files');











%% Create afq structure

if notDefined('out_name')
    out_name = ['afq_', getDateAndTime];
end

disp('Running AFQ_create with the following options...');
fprintf('sub_dirs: %s', sub_dirs{1})
fprintf('output_dir: %s', output_dir)
fprintf('out_name: %s', out_name)
mkdir(input_dir, 'bin')
afq = AFQ_Create('sub_dirs', sub_dirs, 'sub_group', sub_group, ...
                 'outdir', output_dir, 'outname', out_name, ...
                 'params', J.params);  
disp('... end running AFQ_Create')
% disp(afq.params);

% Run control comparison by default
% if ~isfield(params, 'runcontrolcomp');a
%     afq.params.runcontrolcomp = false;
% end





%% RUN AFQ

disp('Running AFQ_run with the following options...');
fprintf('sub_dirs: %s', sub_dirs{1})
disp('This is the afq struct');
afq
afq = AFQ_run(sub_dirs, sub_group, afq);
disp('... end running AFQ_run');

%% Check for empty fiber groups
disp('Checking for empty fiber groups...');
for i = 1:numel(afq.TractProfiles)
    if isempty(afq.TractProfiles(i).nfibers)
        disp(fprintf('Fiber group is empty: %s', afq.TractProfiles(i).name));
    end
end


%% Export the data to csv files (don't use AFQ_exportData)

disp('Exporting data to csv files...');

% We will add the diffusion parameters and the series number to the name
csv_dir = fullfile(output_dir,'csv_files');
mkdir(csv_dir);

% Get the names of each of the datatypes (e.g.,'FA','MD', etc.)
properties = fieldnames(afq.vals);

% Loop over the properties and create a table of values for each property
for ii = 1:numel(properties)

    % Loop over each fibergroup and insert the values into the table
    for i = 1:numel(afq.fgnames)

        % If this is the first time through create a table that we'll
        % concatenate with each following table (t)
        if i == 1
            T = cell2table(num2cell(afq.vals.(properties{ii}){i}'),'variablenames',{regexprep(afq.fgnames{i},' ','_')});
        else
            t = cell2table(num2cell(afq.vals.(properties{ii}){i}'),'variablenames',{regexprep(afq.fgnames{i},' ','_')});

            % Combine the tables
            T = horzcat(T,t);
        end

        % Write out the table to a csv file
        writetable(T,fullfile(csv_dir,['AFQ_' lower(properties{ii}) '.csv']));
    end
end

%% Create the tck files for visualizing the results
%{
We want to see the following things for checking the quality and/or continuing
the processing of the tracts. 
Visualize (all after alignment):
      T1w file
      5tt file (hopefuly based on FS-s aparc+aseg.mgz)
      _fa.mif
      _brainmask.mif
      _wmMask.mif and _wMask_dilated.mif
      _csd file, usually a good idea to check in the first subject to check
          bvec alignment
Visualize final results (first check if the values make sense)
      whole tractogram.tck: how does it fill the WM?
      TRACTS
        - Uncleaned
        - Cleaned
        - ROIs
        
%}

disp('Creating the tck files for visualization and QA...');
% We will add the diffusion parameters and the series number to the name
vis_dir = fullfile(output_dir,'vis_files');
mkdir(vis_dir);

% First of all we will copy all the files present in the dtiInit root to afq so
% that everything is in the same zip
inputParts  = split(input_dir, filesep);
predti = strjoin(inputParts(1:(length(inputParts)-1)), filesep);
copyfile([predti '/*.mat'], output_dir);
copyfile([predti '/*.bv*'], output_dir);
copyfile([predti '/*.nii*'], output_dir);

% Obtain the files
%if isdeployed
    % Convert the ROIs from mat to .nii.gz
    % Read the b0
    img  = niftiRead(fullfile(input_dir, 'bin', 'b0.nii.gz'));
    % Obtain the ROIs in nifti to check if they look ok or not
    rois = dir(fullfile(input_dir, 'ROIs', '*.mat'));
    for df=1:length(rois)
        roiFullPath = fullfile(input_dir, 'ROIs',rois(df).name);
        roi         = dtiReadRoi(roiFullPath);
        coords      = roi.coords;
        % Convert vertex acpc coords to img coords
        imgCoords  = mrAnatXformCoords(img.qto_ijk, coords);
        % Get coords for the unique voxels
        imgCoords = unique(ceil(imgCoords),'rows');
        % Make a 3D image
        roiData = zeros(img.dim);
        roiData(sub2ind(img.dim, imgCoords(:,1), imgCoords(:,2), imgCoords(:,3))) = 1;
        % Change img data
        img.data = roiData;
        img.cal_min = min(roiData(:));
        img.cal_max = max(roiData(:));
        % Write the nifti file
        [~,roiNameWoExt] = fileparts(roiFullPath);
        img.fname = fullfile(fileparts(roiFullPath), [roiNameWoExt,'.nii.gz']); 
        writeFileNifti(img);     
    end
    % Convert the segmented fg-s to tck so that we can see them in mrview
    % In the future we will make them obj so that they can be visualized in FW
    % First create another two MoriSuperFibers out of the clipped and not
    % clipped ones. 
    FGs = dir(fullfile(input_dir, 'fibers', 'Mori*.mat'));
    for nf = 1:length(FGs)
        fgname             = fullfile(input_dir, 'fibers', FGs(nf).name);
        fg                 = fgRead(fgname);
        fgSF               = fg;
        % Change the fiber by the superfiber
        for nsf=1:length(fg)
            if ~isempty(fg(nsf).fibers)
                [SuperFiber] = dtiComputeSuperFiberRepresentation(fg(nsf),[],100);
                fgSF(nsf).fibers= SuperFiber.fibers;
            end
        end
        % Save the clipped ones as well for QA
        [path, fname, fext] = fileparts(fgname);
        sfFgName = fullfile(path,[fname '_SF' fext]);
        dtiWriteFiberGroup(fgSF, sfFgName);
    end
    % Now we will have the same and the newly created ones, that we will
    % create the superFibers.
    FGs = dir(fullfile(input_dir, 'fibers', 'Mori*.mat'));
    for nf = 1:length(FGs)
        fgname             = fullfile(input_dir, 'fibers', FGs(nf).name);
        fg                 = fgRead(fgname);
        saveToMrtrixFolder = true; createObjFiles     = true; 
        AFQ_FgToTck(fg, fgname, saveToMrtrixFolder, createObjFiles)
    end
    % Now we can copy the output to the vis_files, and decide later if copying
    % it to results so that it can be visualized in FW
%{    
else
    % This will be used to download the files when using matlab online:
    % The important part is to make work the isdeployed part. 
    % The code below will be run usually manually, usually for older FW analysis
    % that didn't have the previous code.
    st                    = scitran('stanfordlabs');
    colecName             = '00_VIS';
    analysisLabelContains = 'AllV02: Analysis afq-pipeline-3';
    zipNameContains       = 'AFQ_Output_';
    listOfFilesContain    = {'MoriGroups_clean','_wmMask.mif','_wmMask_dilated.mif', ...
                             '_fa.mif', 'b0.nii.gz','_L.mat','_R.mat'};
    downloadDir           = '/Users/glerma/Downloads/AllV02';
    downFiles             = dr_fwDownloadFileFromZip(st, colecName, zipNameContains, ...
                             'analysisLabelContains', analysisLabelContains, ...
                             'filesContain'         , listOfFilesContain, ...
                             'downloadTo'           , downloadDir, ...
                             'showListSession'      , false);
    % Read template in the same space to write the .mats
    b0Path = downFiles{contains(downFiles,'bin/b0.nii.gz')};
    img = niftiRead(b0Path);
    MoriCleans = {}; nMC = 0;
    for df=1:length(downFiles)
        if contains(downFiles{df}, 'fibers/MoriGroups_clean')
            nMC = nMC + 1;
            MoriCleans{nMC} = downFiles{df};
        end
        if contains(downFiles{df}, 'fibers/MoriGroups')
            fg       = fgRead(downFiles{df});
            saveToMrtrixFolder = true;
            createObjFiles     = false; % we want this in FW, not locally
            AFQ_FgToTck(fg, downFiles{df}, saveToMrtrixFolder, createObjFiles)
        end
        if contains(downFiles{df}, 'ROIs/')
            roi = dtiReadRoi(downFiles{df});
            coords = roi.coords;
            % convert vertex acpc coords to img coords
            imgCoords  = mrAnatXformCoords(img.qto_ijk, coords);
            % get coords for the unique voxels
            imgCoords = unique(ceil(imgCoords),'rows');
            % make a 3D image
            roiData = zeros(img.dim);
            roiData(sub2ind(img.dim, imgCoords(:,1), imgCoords(:,2), imgCoords(:,3))) = 1;
            % change img data
            img.data = roiData;
            img.cal_min = min(roiData(:));
            img.cal_max = max(roiData(:));
            % write the nifti file
            [~,roiNameWoExt] = fileparts(downFiles{df});
            img.fname = fullfile(fileparts(downFiles{df}), [roiNameWoExt,'.nii.gz']); 
            writeFileNifti(img);
        end
    end
    % We need to create the clipped fibers for visualization
    % Copy the same code used inside AFQ to do the same, we will create a
    % new Mori file and the obtain the tck-s as fot the other Mori-groups
    % For this, we need to be sure that all the roi-s have been
    % downloaded previously...
    for nMC=length(MoriCleans)
        fg_clip = fgRead(MoriCleans{nMC});
        dtiDir  = strrep(fileparts(MoriCleans{nMC},'/fibers',''));
        % Remove all fibers that are too long and too far from the core of
        % the group.  This algorithm will constrain the fiber group to
        % something that can be reasonable represented as a 3d gaussian
        for jj = 1:20
            % load ROIs
            [roi1, roi2] = AFQ_LoadROIs(jj,dtiDir);
            fg_clip(jj) = dtiClipFiberGroupToROIs(fg_clean(jj),roi1,roi2);
        end
        % Save the clipped ones as well for QA
        [path, fname, fext] = fileparts(MoriCleans{nMC});
        clippedFgName       = fullfile(path,[fname '_CLIPPED' fext]);
        dtiWriteFiberGroup(fg_clip, clippedFgName);
        % Obtain the superfiber of the cleaned one first
        fgClipSF = fg_clip;
        % Change the fiber by the superfiber
        for nf=1:length(fg)
            SuperFiber = dtiComputeSuperFiberRepresentation(fg_clip(nf),[],100);
            fgClipSF(nf).fibers= SuperFiber.fibers;
        end
        % Save the clipped ones as well for QA
        sfClipFgName = fullfile(path,[fname '_CLIPPED_SF' fext]);
        dtiWriteFiberGroup(fgClipSF, sfClipFgName);
        % In the same place we have the clipped and the clippedSF, create tcks
        AFQ_FgToTck(fg_clip, clippedFgName, saveToMrtrixFolder, createObjFiles)
        AFQ_FgToTck(fgClipSF, sfClipFgName, saveToMrtrixFolder, createObjFiles)
    end
end
%}
%% Create Plots and save out the images
%{
if afq.params.runcontrolcomp

    disp('Running comparison to control population!')

    % Setup the valnames.
    valnames = fieldnames(afq.vals);

    % Remove those values we're not interested in currently
    %TODO: Remove this - generate them all.
    remlist = {'cl','curvature','torsion','volume','WF_map_2DTI','cT1SIR_map_2DTI'};
    for r = 1:numel(remlist)
        ind = cellfind(valnames,remlist{r});
        if ~isempty(ind)
            valnames(ind) = [];
        end
    end

    % Load up saved controls data based on the parameters for the data.
    % Handle the case where ndirs is not 96 or 30.
    a = abs(metadata.ndirs - 96);
    b = abs(metadata.ndirs - 30);

    % If the diffusion values are not exact then we have to warn that there was
    % no exact match! || OR do we just not perform the plotting aginst the
    % control data.
    if metadata.ndirs ~= 96 && metadata.ndirs ~= 30
        warning('Number of diffusion directions does not match control data!');
    end

    if metadata.bvalue ~= 2000 && metadata.bvalue ~= 1000
        warning('B-VALUE does not match control data!');
    end

    if metadata.ndirs == 96 || a < b
        disp('Loading 96-direction control data');
        afq_controls = control_data.afq96;

    elseif metadata.ndirs == 30 || b < a
        disp('Loading 30-direction control data');
        afq_controls = control_data.afq30;
    end

    % This might be where we write out the figures in two seperate directories.
    fig_out_dir = fullfile(output_dir, 'figures');
    if ~exist(fig_out_dir,'dir'); mkdir(fig_out_dir); end

    disp('Running AFQ Plot: Generating figures...');
    try
        if metadata.age_comp && isnumeric(metadata.age)
            disp('Constraining norms based on age!');

            AFQ_PlotPatientMeans(afq, afq_controls, valnames, 21:80, fig_out_dir,'Age', metadata.age_range);
        else
            AFQ_PlotPatientMeans(afq, afq_controls, valnames, 21:80, fig_out_dir);
        end
    catch ME
        disp(ME.message);
    end
end
%}

%% Reproducibility

R = {};
R.date = getDateAndTime;
[~, R.arch] = system('lsb_release -a');
R.software.version = version;
R.software.libs = ver;
R.code = {};

% R.analysis.metadata = metadata;
R.analysis.params   = afq.params;
R.analysis.subject  = sub_dirs;

save(fullfile(output_dir,'Reproduce.mat'), 'R');
savejson('', R, fullfile(output_dir,'Reproduce.json'));


%% END

if isdeployed
    disp('Sending exit(0) signal.');
    exit(0)
else
    disp('Done!');
end

return

% Use compile.sh to compile this file 
% This is the command used to launch it with the MCR
% ./run_AFQ_StandAlone_QMR.sh /data/localhome/glerma/soft/matlab/mcr92/v92   '{\"input_dir\" : \"/data/localhome/glerma/TESTDATA/AFQ/input/dtiInit111_mcr/dti90trilin\", \"output_dir\": \"/data/localhome/glerma/TESTDATA/AFQ/output/withDtiinit111_mrtrix_mcr\",\"params\"    :\"/data/localhome/glerma/TESTDATA/AFQ/input/config_parsed.json\"}'

% After adding LiFE, add this packages too:
%       addpath(genpath('/data/localhome/glerma/soft/encode'));
%       addpath(genpath('/data/localhome/glerma/soft/app-life'));
% So the new mcc command is as follows:
% mcc -m -I /data/localhome/glerma/soft/encode -I /data/localhome/glerma/soft/app-life -I /black/localhome/glerma/soft/spm8 -I /data/localhome/glerma/soft/vistasoft -I /data/localhome/glerma/soft/jsonlab /data/localhome/glerma/soft/afq-pipeline/afq/source/bin/AFQ_StandAlone_QMR.m  

