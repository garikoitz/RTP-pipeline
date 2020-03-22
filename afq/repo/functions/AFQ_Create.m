function afq = AFQ_Create(sub_dirs, params,tractparams,varargin)
% Create AFQ structure and set default parameters
%
%    afq = AFQ_Create(varargin)
%
% Creates an automated fiber quantification (AFQ) structure.  The default
% fields are put in place and filled with default values.  The default
% parameters can also be changed and this will affect later stages of the
% AFQ pipeline.  The arguments to the function are in the form
% 'parameter','value'.  For example, if you call the function with
%
%
% See Also: AFQ_Set AFQ_Get
%
% Examples:
%    afq = AFQ_Create;
%    afq = AFQ_Create('runmode','test');
%
%
% (c) Jason D. Yeatman December 2011
% Garikoitz Lerma Usabiaga 2020



% Parse
p = inputParser;
p.addRequired('sub_dirs', @ischar);
p.addRequired('params', @isstruct);
p.addRequired('tractparams', @istable);

% Parse. Assign result inside each case
p.parse(sub_dirs,params,tractparams,varargin{:});



%% Define the type of structure
afq.type = 'afq version 1';


% Add this for debugging. It won't affect on normal Docker usage because the files won't exist. 
afq.force = false;

% Add the main components. The rest maintain them empty from now. Deleted default parameters below
afq.params = params;
afq.tracts = tractparams;

%% Names of all the fiber groups
afq.fgnames = cellstr(tractparams.label');
%% Names of the ROIs associated with each fiber group
afq.roi1names = cellstr(tractparams.roi1');
afq.roi2names = cellstr(tractparams.roi2');

%% Attach a cell array of subject directories to the afq structure
afq.sub_dirs = {sub_dirs};


% old now
%% Attach a structure of values to the afq structure
vals.fa = {};
vals.md = {};
vals.rd = {};
vals.ad = {};
vals.cl = {};
afq.vals = vals;

%% Attach a cell array of subject ids to the afq structure
afq.sub_ids = {};

%% Attach a vector of subject groups to afq structure
afq.sub_group = 1;

%% Attach a vector of subject id numbers to the afq structure
afq.sub_nums = [];


%% Structured array of scan parameters for data acquisition

% string for vendor. (e.g., 'GE' or 'Siemens')
afq.scanparams.vendor = [];
% string for scanner model (e.g., 'MR750')
afq.scanparams.model = [];
% Field strength if Tesla. Defaults to 3
afq.scanparams.fieldstrength = 3;
% Spatial resolution in mm [x, y, z]. (e.g., [2 2 2])
afq.scanparams.resolution = [];
% List the b-values used for each volume or, at least, the different
% b-values that were used (e.g., 0, 1000, 2000)
afq.scanparams.bvals = [];
% Number of volumes collected at each b-value
afq.scanparams.nvols = [];
% Attatch the full gradient table (e.g., bvecs file)
afq.scanparams.gradtable = [];
% TR for scan
afq.scanparams.TR = [];
% TE for scan
afq.scanparams.TE = [];
% Number of receive coils. (e.g., 8)
afq.scanparams.coils = [];
% If parallel imaging was used (e.g., SENSE) the provide the acceleartion factor
afq.scanparams.acceleration = []; 
% If multiband was used, then provide SMS factor
afq.scanparams.multiband = [];
% If known, provide the diffusion time (or times) associated with each bval
afq.scanparams.diffusiontime = [];

%% Attach the tract profile structure to the afq structure

afq.TractProfiles = AFQ_CreateTractProfile;
% % Add a field to afq.TractProfiles for each tract defined above
% for ii = 1:length(fgNames)
%     afq.TractProfiles.(fgNames{ii}) = AFQ_CreateTractProfile('name',fgNames{ii});
% end

%% Attatch a field for spatial normalization
afq.xform.sn      = [];
afq.xform.invDef  = [];
afq.xform.ants    = [];
afq.xform.antsinv = [];

%% Check which software packages are installed
afq.software.mrvista = check_mrvista;
afq.software.mrtrixVersion = check_mrTrix_Version;
if check_mrTrix_Version ~= 0
    afq.software.mrtrix = check_mrTrix(afq.software.mrtrixVersion);
else
    error('mrTrix not installed or not properly set up')
end



%% Attach a structure pointing to each subjects data files
afq.files.dt6{1} = fullfile(afq.sub_dirs{1},'dt6.mat');
if ~exist(AFQ_get(afq,'dt6 path',1))
    error('%s file does not exist',AFQ_get(afq,'dt6 path',ii))
end
afq.files.images            = struct('name',{},'path',{});
afq.files.fibers.wholebrain = cell(AFQ_get(afq,'num subs'),1);
afq.files.fibers.segmented  = cell(AFQ_get(afq,'num subs'),1);
afq.files.fibers.clean      = cell(AFQ_get(afq,'num subs'),1);


disp('[AFQ_Create] Until here, parameter reading and setup.')
disp('[AFQ_Create] This is afq and afq.files')
afq
afq.files




%% Add files from previous AFQ runs to afq structure

% The name of the segmented fiber group depends on whether we are clipping
% it to the ROIs or not. Or it can be passed in by the user
s = strcmp('segname',varargin) + strcmp('segName',varargin);
if sum(s)>0
    segName = varargin{find(s)+1};
elseif AFQ_get(afq,'clip2rois') == 0
    segName = 'MoriGroups_Cortex.mat';
else
    segName = 'MoriGroups.mat';
end
fibDir       = fullfile(afq.sub_dirs{1},'fibers');
wholebrainFG = fullfile(fibDir,'WholeBrainFG.mat');
% segmentedFG = fullfile(fibDir,segName);
% cleanFG = fullfile(fibDir,[prefix(segName) '_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']);
if exist(wholebrainFG,'file')
    afq.files.fibers.wholebrain{1} = wholebrainFG;
end
% Save the name  of the segmented fiber group
afq.files.fibers.segName = segName;

%% Allow previous analyses to be overwritten
afq.overwrite.fibers.wholebrain = zeros(AFQ_get(afq,'num subs'),1);
afq.overwrite.fibers.segmented = zeros(AFQ_get(afq,'num subs'),1);
afq.overwrite.fibers.clean = zeros(AFQ_get(afq,'num subs'),1);
afq.overwrite.vals = zeros(AFQ_get(afq,'num subs'),1);

%% If desired compute constrained spherical deconvolution with mrtrix
% If we want to use mrtrix for tractography that we will compute CSD right
% here


% For some reason it is not going inside this loop, at least when testing in black the docker container
%if AFQ_get(afq,'use mrtrix')
ii = 1;
mrtrixdir = fullfile(afq.sub_dirs{ii},'mrtrix');
if ~exist(mrtrixdir,'dir'),mkdir(mrtrixdir);end
% Get the lmax from the afq structure
% lmax = AFQ_get(afq,'lmax');
% Obtain the lmax from the bvalues, no need to set up, calculate
if afq.params.track.mrtrix_autolmax
    lmax = AFQ_calculateLmaxFrombvals(AFQ_get(afq, 'dt6path',ii));
else
    lmax = afq.params.track.mrtrix_lmax;
end
% Beware, in the latest versions of mrtrix we are not using lmax,
% their scripts take care of it. 
% So, there is manual, there is autoMax (calculated above) and
% automrtrix, calculated by them Right now, the code will just do
% automrtrix... Change it and leave it explained. 
dt6path =AFQ_get(afq, 'dt6path',ii) ;
fprintf('This is the dt6path going to AFQ_mrtrixInit: %s',dt6path);
files = AFQ_mrtrixInit(dt6path, ...
                       lmax,...
                       mrtrixdir,...
                       afq.software.mrtrixVersion, ...
                       afq.params.track.multishell, ... % true/false
                       afq.params.track.tool, ... % 'fsl', 'freesurfer'
                       afq.params.track.faMaskThresh, ...
                       afq.force);
% In order to not modify much the previous code, I created new
% files types. 
% In mrTrix2 and mrTrix3 not-multishell, files.wm was the wm mask,
% so I changed the name to files.wmMask.
% In multishell, in files.tt5 you have the wm, gm, csf masks in one
% file. We create it only if it is multishell, but wmMask is always
% created because we will need it downstream in tractography.
% files.csd is created  only in  ~multishell and passed here to
% tractography, but in the case of msmt 3 different files are
% created, one for each tissue type. We only pass the csd of the 
% wm = wmMask for tractography, wmMask as seed_image
% and tt5 for -act (instead of -mask)
afq.files.mrtrix.wm{ii}         = files.wmMask;
afq.files.mrtrix.wm_dilated{ii} = files.wmMask_dilated;
afq.files.mrtrix.tt5{ii}        = files.tt5;
afq.files.mrtrix.gmwmi{ii}      = files.gmwmi;
% Now we are always providing the wmCsd to the tractography
afq.files.mrtrix.csd{ii}        = files.wmCsd;
% if afq.params.track.multishell; afq.files.mrtrix.csd{ii} = files.wmCsd;
% else; afq.files.mrtrix.csd{ii}   = files.csd;
% end

% This is new, here now we will create the files for the /bin
% folder that we did not create in AFQ_dtiInit.m
bindir = fullfile(afq.sub_dirs{ii},'bin');
% These are the files to be created:
binfiles.b0        = fullfile(bindir,'b0.nii.gz');
binfiles.brainmask = fullfile(bindir,'brainMask.nii.gz');
binfiles.wmMask    = fullfile(bindir,'wmMask.nii.gz');
binfiles.tensors   = fullfile(bindir,'tensors.nii.gz');
binfiles.fa        = fullfile(bindir,'fa.nii.gz');

% use a series of AFQ_mrtrix_convert-s for this
if exist(binfiles.b0,'file') && afq.force == false
    fprintf('[AFQ_Create] %s exists and will not overwritten\n',binfiles.b0);
else; AFQ_mrtrix_mrconvert(files.b0, binfiles.b0,0,0,afq.software.mrtrixVersion);end
% The b0 coming from mrtrix has multiple volumes, and I think
% mrDiffusion is expecting a single volume, obtain the mean here
A       = niftiRead(binfiles.b0); % Reading nifti created by mrtrix
A.data  = mean(A.data,4);
A.dim   = size(A.data);
A.ndim  = length(A.dim);
A.pixdim= A.pixdim(1:3);
niftiWrite(A);
% Continue with the rest of conversions
% disp('Here the conversion of mif files to /bin/*.nii.gz-s is done');
if exist(binfiles.brainmask,'file') && afq.force == false
    fprintf('[AFQ_Create] %s exists and will not be overwritten\n',binfiles.b0);
else; AFQ_mrtrix_mrconvert(files.brainmask, binfiles.brainmask,0,0,afq.software.mrtrixVersion); end

if exist(binfiles.wmMask,'file') && afq.force == false
    fprintf('[AFQ_Create] %s exists and will not be overwritten\n',binfiles.wmMask);
else; AFQ_mrtrix_mrconvert(files.wmMask, binfiles.wmMask,0,0,afq.software.mrtrixVersion); end

if exist(binfiles.tensors,'file') && afq.force == false
    fprintf('[AFQ_Create] %s exists and will not be overwritten\n',binfiles.tensors);
else; AFQ_mrtrix_mrconvert(files.dt, binfiles.tensors,0,0,afq.software.mrtrixVersion); end

if exist(binfiles.fa,'file') && afq.force == false
    fprintf('[AFQ_Create] %s exists and will not be overwritten\n',binfiles.fa);
else; AFQ_mrtrix_mrconvert(files.fa, binfiles.fa,0,0,afq.software.mrtrixVersion);  end

% In order to make the rest of the flow work well, we will modify the
% tensor file to be the same mrDiffusion is expecting
B       = niftiRead(binfiles.tensors); % Reading nifti created by mrtrix
sz      = size(B.data);
if sz(4)~=1
    B.data  = reshape(B.data,[sz(1:3),1,sz(4)]);
    % This is horrible. Mrtrix and mrDiffusion use the same format dxx,dyy,dzz...
    % in order to make the workflow work I need to convert it when
    % writing and when reading in order to maintain existing functions
    B.data = B.data(:,:,:,1,[1 4 2 5 6 3]);
    B.dim   = size(B.data);
    B.ndim  = length(B.dim);
    B.pixdim= [B.pixdim, 1];
    niftiWrite(B);
end
         

%% Set the current subject field to subject 1
afq.currentsub = 1;

%% Add a field for meta data (eg. age, sex etc.)
afq.metadata = [];

return
