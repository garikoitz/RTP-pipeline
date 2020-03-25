function RTP(jsonargs)
%
% RTP_StandAlone_QMR(jsonargs)
%
% INPUT ARGUMENTS: Note: must be a json string.
%
%       input_dir:  Directory containing dt6.mat
%       out_name:   Name for the resulting afq file
%       output_dir: Location to save the results
%       params:     Key-Value pair of params for RTP_Create
%{
% EXAMPLE USAGE:
jsonargs = "/black/localhome/glerma/TESTDATA/FS/paramsBLACK.json";
jsonargs = "/data/localhome/glerma/soft/RTP-pipeline/example_output_parsed.json";
RTP(jsonargs);
%}
%
%#ok<*AGROW>


%% Begin

disp('[RTP] Starting RTP...');
datestamp = strrep(char(datetime(datestr(now),'TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ssz')),':','_');
disp(['[RTP] datestamp:', datestamp]);

% Initialize args
input_dir  = [];
output_dir = [];
out_name   = [];
params     = [];
anat_dir   = [];
bval_dir   = [];
bvec_dir   = [];
fs_dir     = [];
nifti_dir  = [];
tractparams_dir = [];

%% Handle jsonargs
disp('[RTP] This is the json string to be read by loadjson:')
disp(jsonargs)


% Read the file in the variable P
if exist('jsonargs', 'var') && ~isempty(jsonargs)
    P = jsondecode(fileread(jsonargs));
end
P
%% Configure inputs and defaults, get rid of P:
input_dir  = P.input_dir;
output_dir = P.output_dir;
params     = P.params;
out_name   = ['rtp_', datestamp];
params     = P.params;
anat_dir   = P.anat_dir;
bval_dir   = P.bval_dir;
bvec_dir   = P.bvec_dir;
fs_dir     = P.fs_dir;
nifti_dir  = P.nifti_dir;
tractparams_dir = P.tractparams_dir;



% Copy files from input dir (probably read only) to the output dir
rtp_dir    = fullfile(output_dir,'RTP');
if exist(rtp_dir)
	error('[RTP] rtp_dir exists in %s', rtp_dir)
else
	mkdir(rtp_dir);
end
% We need these files in root dir (rtp_dir) to start working
t1_file    = fullfile(rtp_dir, 't1.nii.gz');
bvec_file  = fullfile(rtp_dir, 'dwi.bvecs');
bval_file  = fullfile(rtp_dir, 'dwi.bvals');
dwi_file   = fullfile(rtp_dir, 'dwi.nii.gz');
fs_file    = fullfile(rtp_dir, 'fs.zip');
%% Copy input files to output/RTP/
copyfile(fullfile(anat_dir ,'t1.nii.gz' ), t1_file)  ;
copyfile(fullfile(bvec_dir ,'dwi.bvecs' ), bvec_file);
copyfile(fullfile(bval_dir ,'dwi.bvals' ), bval_file);
copyfile(fullfile(nifti_dir,'dwi.nii.gz'), dwi_file) ;
copyfile(fullfile(fs_dir   ,'fs.zip'    ), fs_file)  ;
% Check if copied properly
if ~exist(t1_file);  error('%s file is not there', t1_file);end
if ~exist(bvec_file);error('%s file is not there', bvec_file);end
if ~exist(bval_file);error('%s file is not there', bval_file);end
if ~exist(dwi_file); error('%s file is not there', dwi_file);end
if ~exist(fs_file);  error('%s file is not there', fs_file);end

% Unzip the output from FS
unzip(fs_file, rtp_dir)
% Check if the expected folders are here, and the minumun files so that the rest does not fail
if ~exist(fullfile(rtp_dir,'fs'),'dir');error('fs folder does not exist, check the fs.zip file');end
if ~exist(fullfile(rtp_dir,'fs','ROIs'),'dir');error('fs/ROIs folder does not exist, check the fs.zip file');end
if ~exist(fullfile(rtp_dir,'fs','aparc+aseg.nii.gz'),'file');error('fs/ROIs/aparc+aseg.nii.gz file does not exist, check the fs.zip file');end
if ~exist(fullfile(rtp_dir,'fs','brainmask.nii.gz'),'file');error('fs/ROIs/brainmask.nii.gz file does not exist, check the fs.zip file');end


% RTP dir will be zipped at the end, the whole dir
% output_dir will be the output in FW, it will have the RTP dir and the files we want to be accesible for reading in FW output
% Only one subject, but until changed this is a cell array
sub_dirs{1} = rtp_dir;


%% Initialize diffusion parameters. This was for dtiInit and dt6Creation. Not really required. TODO: remove this
dwParams            = dtiInitParams;
dwParams.outDir     = output_dir;
dwParams.bvecsFile  = bvec_file;
dwParams.bvalsFile  = bval_file;

%% Validate that the bval values are normalized
% Determine shell
bvals        = dlmread(bval_file);
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
dlmwrite(bval_file, roundedBval, 'delimiter',' ');



%% Run (the equivalent ) dtiInit
% From 3.0.5 onwards I am forking dtiInit and giving it less
% functionalities. I will remove the tensor fitting and do it with mrTrix,
% it is much faster and the rest of the code relies in it anyways. 
% First iteration, stop dtiInit doing it and later on mrtrixInit will paste
% the fa to the bin folder. 
% 3.1.1 it was doing nothing already
% 3.1.2 I am removing the call altogether to make this thing simpler
%        AFQ_dtiInit(J.dwi_file, J.t1_file, dwParams);
% 3.2.0: remove afq-browser and dtiInit altogether
% 4.0.0 now, that we already removed dtiInit and afq-browser.
% 4.2.0 reads csv and creates tracts table, now any tract can be done
% 4.2.1 ROIs: dilates and concatenates them

% Here dtiInit was called, assign variables
dwRawFileName = dwi_file;
t1FileName    = t1_file;

% I. Load the diffusion data, set up parameters and directories structure
% Load the difusion data
disp('Loading preprocessed (use rtp-preproc or already preprocessed) data...');
dwRaw = niftiRead(dwRawFileName);
dwParams.dwOutMm = dwRaw.pixdim(1:3);

% Initialize the structure containing all directory info and file names
% dwDir      = dtiInitDir(dwRawFileName,dwParams);
dwDir.outSuffix        = ''
dwDir.mrDiffusionDir   = mrdRootPath();
dwDir.dataDir          = rtp_dir;
dwDir.inBaseName       = 'dwi'
dwDir.subjectDir       = dwDir.dataDir;
dwDir.mnB0Name         = fullfile(dwDir.dataDir, 'dwi_b0.nii.gz');
dwDir.outBaseName      = ''
dwDir.outBaseDir       = dwDir.dataDir
dwDir.inBaseDir        = fullfile(dwDir.dataDir,'dwi');
dwDir.bvalsFile        = fullfile(dwDir.dataDir,'dwi.bvals');
dwDir.bvecsFile        = fullfile(dwDir.dataDir,'dwi.bvecs');
dwDir.ecFile           = fullfile(dwDir.dataDir,'dwi_ecXform.mat');
dwDir.acpcFile         = fullfile(dwDir.dataDir,'dwi_acpcXform.mat');
dwDir.alignedBvecsFile = dwDir.bvecsFile;
dwDir.alignedBvalsFile = dwDir.bvalsFile;
dwDir.dwAlignedRawFile = dwi_file;

outBaseDir = dwDir.outBaseDir;
fprintf('Dims = [%d %d %d %d] \nData Dir = %s \n', size(dwRaw.data), dwDir.dataDir);
fprintf('Output Dir = %s \n', dwDir.subjectDir);

dwParams.dt6BaseName = fullfile(rtp_dir);
dt6FileName          = fullfile(rtp_dir, 'dt6.mat');
binDirName           = fullfile(rtp_dir, 'bin');
if(~exist(binDirName,'dir')) ;mkdir(binDirName);end
if(~exist('adcUnits','var')); adcUnits = ''; end


params.input_dir    = rtp_dir;
params.output_dir   = output_dir;
params.buildDate    = datestamp; 
l = license('inuse');
params.buildId      = sprintf('%s on Matlab R%s (%s)',l(1).user,version('-release'),computer);
if(ischar(dwRawFileName)); [dataDir,rawDataFileName] = fileparts(dwRawFileName);  % dwRawAligned);
else                       [dataDir,rawDataFileName] = fileparts(dwRawFileName.fname);end  % dwRawAligned.fname); end
params.rawDataDir   = dataDir;
params.rawDataFile  = rawDataFileName;
params.subDir       = dataDir;
params.fs_dir       = fullfile(rtp_dir,'fs');
params.roi_dir      = fullfile(rtp_dir,'ROIs');
% These need to be relative to output_dir...
files.b0            = fullfile('RTP','bin','b0.nii.gz');
files.brainMask     = fullfile('RTP','bin','brainMask.nii.gz');
files.wmMask        = fullfile('RTP','bin','wmMask.nii.gz');
files.tensors       = fullfile('RTP','bin','tensors.nii.gz');
files.fa            = fullfile('RTP','bin','fa.nii.gz');

save(dt6FileName,'adcUnits','params','files');
dtiInitDt6Files(dt6FileName,dwDir,t1FileName);


%% CHECK IF INPUT IS RAS  (THIS WAS IN RTP CONTAINER)
% New in 3.1.2: check the file is RAS If not, convert to RAS Be careful, bvecs
% needs to be changed as well accordingly Instead of changing bvecs manually, we
% will let mrtrix take care of it Convert files to mif, add the bvecs and bvals
% as part of the file, do the conversion and then output the fsl type bvecs
% again. This we can maintain coherence in the whole process. Therefore, if
% there is a - in the strides or the order is not 123, we need to convert (with
% freesurfer it is easier because it gives you RAS or PIR or LAS or whatever but
% it is not installed in the Docker container)

% Read the input file names and convert them
basedir = rtp_dir;

fprintf('This is basedir for RAS checks: %s\n',   basedir)

J       = load(fullfile(rtp_dir,'dt6.mat'));
disp('These are the contents of dt6.mat')
J.params.fs_dir = fullfile(rtp_dir,'fs');
J.params.roi_dir = fullfile(rtp_dir,'fs','ROIs');
J.params
J.files

% DWI file
[p,f,e]= fileparts(J.files.alignedDwRaw);
if ~strcmp(p,basedir); J.files.alignedDwRaw = fullfile(basedir,[f e]); end

% BVEC file
[p,f,e]= fileparts(J.files.alignedDwBvecs);
if ~strcmp(p,basedir); J.files.alignedDwBvecs = fullfile(basedir,[f e]); end

% BVAL file
[p,f,e]= fileparts(J.files.alignedDwBvals);
if ~strcmp(p,basedir); J.files.alignedDwBvals = fullfile(basedir,[f e]); end


J.files.t1path    = fullfile(J.params.output_dir,J.files.t1);
fprintf('This is the absolute path to the t1: %s\n', J.files.t1path)
if exist(J.files.t1path,'file')
    fprintf('T1 file %s Exists. \n', J.files.t1path)
else
    error('Cannot find %s', J.files.t1path)
end

% Check it in these files
checkfiles = {'alignedDwRaw','t1path'}; %,'aparcaseg'};
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


% Read the csv with the tract, so that we can pass it to AFQ_Create
% TODO: shorten this script creating independent functions, for example, RAS check, or ROIs check

csv = dir(fullfile(tractparams_dir,'*.csv'));
if length(csv) > 1; 
    error('There are more than one csv files for tract params');
else
    if isfile(fullfile(csv.folder, csv.name));
        A = readtable(fullfile(csv.folder, csv.name), ...
					 'FileType', 'text', ...
					 'Delimiter','comma', ...
					 'ReadVariableNames',true,...
					 'TextType', 'string');
    else
        error('Cannot read %s or is not a file',fullfile(csv.folder, csv.name))
    end
end
% Check that all ROIs are available in the fs/ROIs folder, if not, throw error. 
% Create list of all ROIs
roi3names=[];
for nr=1:length(A.roi3)
	if ismissing(A.roi3(nr))
		roi3names = [roi3names;""];
	else
		roi3names = [roi3names;strcat(A.roi3(nr),A.extroi3(nr))];
	end
end

checkTheseRois = [strcat(A.roi1,A.extroi1);strcat(A.roi2,A.extroi2);roi3names];
% If there is an AND, this means that there are two ROIs 
% Create a new list with the individual ROIs to create, after checking the individuals are there
createROInew = [];
createROI1   = [];
createROI2   = [];
for nc=1:length(checkTheseRois)
	rname = checkTheseRois(nc);
	if strcmp(rname,"NO.nii.gz")
		% do nothing
	elseif contains(rname,'_AND_')
		% Add the ROI to be created
		createROInew = [createROInew; rname];	
		% Check if the individuals exist
		rois12 = strsplit(rname,'_AND_');

		% ROI1				
		rpath  = fullfile(J.params.roi_dir,strcat(rois12{1},".nii.gz"));
		if ~isfile(rpath);error('ROI %s is required and it is not in the ROIs folder',rpath);end
		createROI1 = [createROI1; strcat(rois12{1},".nii.gz")];	
	
		% ROI2
		rpath  = fullfile(J.params.roi_dir,rois12{2});
		if ~isfile(rpath);error('ROI %s is required and it is not in the ROIs folder',rpath);end
		createROI2 = [createROI2; string(rois12{2})];	
	else
		rpath  = fullfile(J.params.roi_dir,rname);
		if ~isfile(rpath);error('ROI %s is required and it is not in the ROIs folder',rpath);end
	end
end
% Create the ROIs by concatenating. Use Matlab for now
for nt=1:length(createROInew)
	nroi = fullfile(J.params.roi_dir,createROInew(nt));
	roi1 = fullfile(J.params.roi_dir,createROI1(nt));
	roi2 = fullfile(J.params.roi_dir,createROI2(nt));
	% Read the existing ROIs
	R1   = niftiRead(char(roi1));
	R2   = niftiRead(char(roi2));
	% Create the new file and concatenate the data
	nR   = R1;
	nR.fname = char(nroi);
	nR.data  = uint8(R1.data | R2.data);
	niftiWrite(nR);  
end

% Dilate those ROIs that require it using mrtrix's tool maskfilter
% maskfilter input filter(dilate) output 
for nl=1:height(A)	
	% Check line by line and create the dilated ROIs.  
	ts = A(nl,:);
	if ts.dilroi1>0
		inroi  = char(fullfile(J.params.roi_dir, strcat(ts.roi1,'.nii.gz')));
		outroi = char(fullfile(J.params.roi_dir, strcat(ts.roi1,'_dil-',num2str(ts.dilroi1),'.nii.gz')));
		cmd    = ['maskfilter -force -npass ' num2str(ts.dilroi1)  ' ' inroi  ' dilate - | '...
				  'mrthreshold -force -abs 0.5 - ' outroi];
		cmdr   = AFQ_mrtrix_cmd(cmd);
		if cmdr ~= 0
			error('[RTP] ROI could not be created, this was the command: %s', cmd)
		end
	end
	if ts.dilroi2>0
		inroi  = char(fullfile(J.params.roi_dir, strcat(ts.roi2,'.nii.gz')));
		outroi = char(fullfile(J.params.roi_dir, strcat(ts.roi2,'_dil-',num2str(ts.dilroi2),'.nii.gz')));
		cmd    = ['maskfilter -force -npass ' num2str(ts.dilroi2)  ' ' inroi  ' dilate - | '...
				  'mrthreshold -force -abs 0.5 - ' outroi];
		cmdr   = AFQ_mrtrix_cmd(cmd);
		if cmdr ~= 0
			error('[RTP] ROI could not be created, this was the command: %s', cmd)
		end
	end
	if ts.dilroi3>0 && ~strcmp(ts.roi3,"NO")
		inroi  = char(fullfile(J.params.roi_dir, strcat(ts.roi3,'.nii.gz')));
		outroi = char(fullfile(J.params.roi_dir, strcat(ts.roi3,'_dil-',num2str(ts.dilroi3),'.nii.gz')));
		cmd    = ['maskfilter -force -npass ' num2str(ts.dilroi3)  ' ' inroi  ' dilate - | '...
				  'mrthreshold -force -abs 0.5 - ' outroi];
		cmdr   = AFQ_mrtrix_cmd(cmd);
		if cmdr ~= 0
			error('[RTP] ROI could not be created, this was the command: %s', cmd)
		end
	end
end
% FINISHED ROI CHECK AND CREATION





%% Create afq structure

disp('Running AFQ_create with the following options...');
fprintf('sub_dirs: %s', sub_dirs{1})
fprintf('output_dir: %s', output_dir)
fprintf('out_name: %s', out_name)
if ~exist(fullfile(input_dir,'bin'));mkdir(input_dir, 'bin');end
% Add deleted variables to afq, the ones we want fixed
J.params.clip2rois       = true;
J.params.track.multishell= false;
if length(paramsShells) > 1;J.params.track.multishell=true;end
J.params.track.tool      = 'freesurfer';
J.params.track.algorithm = 'mrtrix';
J.params.computeCSD      = true;
afq = AFQ_Create(sub_dirs{1}, J.params, A);  



%{
afq.tracts = A;
% Update afq values with the ones coming from the tract. 
% TODO: do this inside AFQ_Create
afq.fgnames = cellstr(afq.tracts.label');
afq.roi1names = cellstr(strcat(afq.tracts.roi1,afq.tracts.extroi1)');
afq.roi2names = cellstr(strcat(afq.tracts.roi2,afq.tracts.extroi2)');
%}

disp('[RTP] ... end running AFQ_Create')

%% RUN RTP

disp('[RTP] Running AFQ_run with the following options...');
fprintf('sub_dirs: %s', sub_dirs{1})
disp('[RTP] This is the afq struct going to AFQ_run');
afq
afq = AFQ_run(sub_dirs, 1, afq);
disp('      ... end running AFQ_run');

%% Check for empty fiber groups
disp('[RTP] Checking for empty fiber groups...');
for i = 1:numel(afq.TractProfiles)
    if isempty(afq.TractProfiles(i).nfibers)
        disp(fprintf('Fiber group is empty: %s', afq.TractProfiles(i).name));
    end
end

% Write the afq file 
% Save each iteration of afq run if an output directory was defined
outname = fullfile(rtp_dir,['afq_' datestamp]);
save(outname,'afq');


%% Export the data to csv files (don't use AFQ_exportData)
disp('[RTP] Exporting data to csv files...');
% We will add the diffusion parameters and the series number to the name
csv_dir = fullfile(rtp_dir,'csv_files');
if ~exist(csv_dir,'dir');mkdir(csv_dir);end

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
        writetable(T,fullfile(csv_dir,['RTP_' lower(properties{ii}) '.csv']));
    end
end



disp('[RTP] Exporting tck files for QA...');
% We will add the diffusion parameters and the series number to the name
tck_dir    = fullfile(rtp_dir,'tck_files');
if ~exist(tck_dir,'dir');mkdir(tck_dir);end
% Copy only the tck files we want to visualize, so that we do not need to unzip the big RTP file
mrtrixdir  = fullfile(J.params.subDir,'mrtrix');
clean_tcks = dir(fullfile(mrtrixdir,'*clean.tck'));
for nt=1:length(clean_tcks)
	srctckname = fullfile(mrtrixdir, clean_tcks(nt).name);
	dsttckname = fullfile(tck_dir  , clean_tcks(nt).name);
	dstplyname = fullfile(tck_dir  , strrep(clean_tcks(nt).name,'tck','ply'));
	copyfile(srctckname,dsttckname)
	% Use the same step to create the ply file in the vis folder. We wil convert them to obj with python
	cmd  = sprintf('tckconvert -dec %s %s', dsttckname, dstplyname);
	rcmd = AFQ_mrtrix_cmd(cmd);
end

%TODO: add a flag, saveIntermediateFilesForQC 
% It will show the files showns below in the results folder

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


% Obtain the files
%if isdeployed
	% disp('Creating the obj files for visualization and QA in FW...');
	% % We will add the diffusion parameters and the series number to the name
	% vis_dir = fullfile(output_dir,'vis_files');
	% mkdir(vis_dir);
    % Convert the ROIs from mat to .nii.gz
    % Not now, now the ROIs come from FS, there is no .mat ROIs anymore
    % Read the b0
    % img  = niftiRead(fullfile(input_dir, 'bin', 'b0.nii.gz'));
    % Convert the segmented fg-s to tck so that we can see them in mrview
    % In the future we will make them obj so that they can be visualized in FW
    % First create another two MoriSuperFibers out of the clipped and not
    % clipped ones. 
    %FGs = dir(fullfile(input_dir, 'fibers', 'Mori*.mat'));
    %for nf = 1:length(FGs)
    %    fgname             = fullfile(input_dir, 'fibers', FGs(nf).name);
    %    fg                 = fgRead(fgname);
    %    fgSF               = fg;
    %    % Change the fiber by the superfiber
    %    for nsf=1:length(fg)
    %        if ~isempty(fg(nsf).fibers)
    %            [SuperFiber] = dtiComputeSuperFiberRepresentation(fg(nsf),[],100);
    %            fgSF(nsf).fibers= SuperFiber.fibers;
    %        end
    %    end
    %    % Save the clipped ones as well for QA
    %    [path, fname, fext] = fileparts(fgname);
    %    sfFgName = fullfile(path,[fname '_SF' fext]);
    %    dtiWriteFiberGroup(fgSF, sfFgName);
    %end
    %% Now we will have the same and the newly created ones, that we will
    %% create the superFibers.
    %FGs = dir(fullfile(input_dir, 'fibers', 'Mori*.mat'));
    %for nf = 1:length(FGs)
    %    fgname             = fullfile(input_dir, 'fibers', FGs(nf).name);
    %    fg                 = fgRead(fgname);
    %    saveToMrtrixFolder = true; createObjFiles     = true; 
    %    AFQ_FgToTck(fg, fgname, saveToMrtrixFolder, createObjFiles)
    %end
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

R = {};
R.date = datestamp;
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

