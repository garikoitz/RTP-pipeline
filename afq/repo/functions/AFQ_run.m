function [afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq)
% Run AFQ analysis on a set of subjects to generate Tract Profiles of white
% matter properties.
%
% [afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs,
% sub_group, [afq])
%
% AFQ_run is the main function to run the AFQ analysis pipeline.  Each AFQ
% function is an independent module that can be run on its own.  However
% when AFQ_run is used to analyze data all the proceeding analyses are
% organized into the afq data structure. The AFQ analysis pipeline is
% described in Yeatman J.D., Dougherty R.F., Myall N.J., Wandell B.A.,
% Feldman H.M. (2012). Tract Profiles of White Matter Properties:
% Automating Fiber-Tract Quantification. PLoS One.
%
% Input arguments:
%  sub_dirs  = 1 x N cell array where N is the number of subjects in the
%              study. Each cell should contain the full path to a subjects
%              data directory where there dt6.mat file is.
%
%  sub_group = Binary vector defining each subject's group. 0 for control
%              and 1 for patient.
%
%  afq       = This is a structure that sets up all the parameters for the
%              analysis.  If it is blank AFQ_run will use the default
%              parameters.  See AFQ_Create.
%
% Outputs: 
% afq          = afq structure containing all the results
%
% patient_data = A 1X20 structured array of tract diffusion profiles where
%                data for each tract is in a cell of the structure (eg.
%                patient_data(1) is data for the left thalamic radiation).
%                Each diffusion properties is stored as a different field
%                (eg. patient_data(1).FA is a matrix of FA profiles for the
%                left thalamic radiation). Within the data matrix each
%                subject is a row and each location is a column.  This
%                output variable contains all the data for the patients
%                defined by sub_group ==1.
%
% control_data = The same structure as for patient_data but this contains
%                data for the control subjects defined by sub_group==0.
%
% norms        = Means and standard deviations for each tract diffusion
%                profile calculated based on the control_data.
%
% abn          = A 1 x N vector where N is the number of patients.
%                Each patient that is abnormal on at least one tract is
%                marked with a 1 and each subject that is normal on every
%                tract is marked with a 0. The criteria for abnormal is
%                defined in afq.params.cutoff.  See AFQ create
%
% abnTracts    = An M by N matrix where M is the number of subjects and N
%                is the number of tracts. Each row is a subject and each 
%                column is a tract.  1 means that tract was abnormal for
%                that subject and 0 means it was normal.
%
%  Web resources
%    http://white.stanford.edu/newlm/index.php/AFQ
%
%  Example:
%   
%   % Get the path to the AFQ directories
%   [AFQbase AFQdata] = AFQ_directories;
%   % Create a cell array where each cell is the path to a data directory
%   sub_dirs = {[AFQdata '/patient_01/dti30'], [AFQdata '/patient_02/dti30']...
%   [AFQdata '/patient_03/dti30'], [AFQdata '/control_01/dti30']...
%   [AFQdata '/control_02/dti30'], [AFQdata '/control_03/dti30']};
%   % Create a vector of 0s and 1s defining who is a patient and a control
%   sub_group = [1, 1, 1, 0, 0, 0]; 
%   % Run AFQ in test mode to save time. No inputs are needed to run AFQ 
%   % with the default settings. AFQ_Create builds the afq structure. This
%   % will also be done automatically by AFQ_run if the user does not wish 
%   % to modify any parameters
%   afq = AFQ_Create('run_mode','test', 'sub_dirs', sub_dirs, 'sub_group', sub_group); 
%   [afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq)
%
% Copyright Stanford Vista Team, 2011. Written by Jason D. Yeatman,
% Brian A. Wandell and Robert F. Dougherty

%% Check Inputs
if notDefined('sub_dirs') && exist('afq','var') && ~isempty(afq)
    sub_dirs = AFQ_get(afq,'sub_dirs');
elseif notDefined('sub_dirs')
    error('No subject directories');
end
if ~iscell(sub_dirs), sub_dirs = cellstr(sub_dirs); end
if notDefined('sub_group') && exist('afq','var') && ~isempty(afq)
    sub_group = AFQ_get(afq,'sub_group');
elseif notDefined('sub_group')
    error('Must define subject group');
end
if length(sub_group) ~= size(sub_dirs,1) && length(sub_group) ~= size(sub_dirs,2)
    error('Mis-match between subject group description and subject data directories');
end
if ~exist('afq','var') || isempty(afq)
    %if no parameters are defined use the defualts
    afq = AFQ_Create('sub_dirs',sub_dirs,'sub_group',sub_group); 
end
if isempty(afq.sub_group)
    afq = AFQ_set(afq,'sub_group',sub_group);
end
if isempty(afq.sub_dirs)
    afq = AFQ_set(afq,'sub_dirs',sub_dirs);
end
% Check which subjects should be run
runsubs = AFQ_get(afq,'run subjects');
% Define the name of the segmented fiber group
segName = AFQ_get(afq,'segfilename');

% If ANTS is installed on the system then precompute spatial normalization
% with ANTS and save to the afq structure
% BW suggests: move ants outputs to ANTs subdir

tic
if AFQ_get(afq, 'use ANTS')
    disp('ANTs normalization, it can take hours')
    afq = AFQ_ComputeSpatialNormalization(afq);
end
toc

%%  Loop over every subject
for ii = runsubs
    % Define the current subject to process
    afq = AFQ_set(afq,'current subject',ii);
    
    %% Preprocess Data
    
    if ~exist(fullfile(sub_dirs{ii},'dt6.mat'),'file')
        error('AFQ preprocessing has not yet been implemented')
    end
    % Load Dt6 File
    dtFile = fullfile(sub_dirs{ii},'dt6.mat');
    dt     = dtiLoadDt6(dtFile);
   
	baseDir = sub_dirs{ii};
    fibDir  = fullfile(baseDir,'fibers');
    mrtrixDir = fullfile(baseDir,'mrtrix');


 
    % If ANTS was used to compute a spatial normalization then load it for
    % this subject
    % EDIT GLU: no, delete this line
    % antsInvWarp = AFQ_get(afq,'ants inverse warp',ii);
    
    %% Perform Whole Brain Streamlines Tractography
    % If required, it will be done in RTP_TractsGet()
    %{
    %Check if there is a fibers directory, otherwise make one.
    fibDir = fullfile(sub_dirs{ii},'fibers');
    if ~exist(fibDir,'dir')
        mkdir(fibDir);
    end 
    % Check if wholebrain tractography should be done
    if AFQ_get(afq, 'do tracking', ii) == 1
        % Perform whole-brain tractography
        fprintf('\nPerforming whole-brain tractograpy for subject %s\n',sub_dirs{ii});
        
        % TODO GLU: here we should add an option to clean data with LiFE or
        % not. Now it will do it by default. fg will be the fibers after
        % LiFE have been run and the 0 weighted fibers have been removed
        fg = AFQ_WholebrainTractography(dt, afq.params.run_mode, afq);
        % Save fiber group to fibers directory
        dtiWriteFiberGroup(fg,fullfile(fibDir,'WholeBrainFG.mat'));
        % Set the path to the fibers in the afq structure
        afq = AFQ_set(afq,'wholebrain fg path','subnum',ii,fullfile(fibDir,'WholeBrainFG.mat'));
        % Wholebrain fiber group is already in memory and does not need to
        % be loaded
        loadWholebrain = 0;
    else
        fprintf('\nWhole-brain tractography was already done for subject %s',sub_dirs{ii});
        % Wholebrain fiber group needs to be loaded
        loadWholebrain = 1;
    end
    %}
  
    %% Segment 20 Fiber Groups
    
    % Check if fiber group segmentation was already done
    % if AFQ_get(afq, 'do segmentation',ii) == 1
        % Load the wholebrain fiber group if necessary
        % if loadWholebrain == 1
        %     fg = AFQ_get(afq,'wholebrain fiber group',ii);
        % end
        % Segment fiber group
        % TODO: delete this line
        % fg_classified = RTP_TractsGet(dtFile, fg, [], [],[], antsInvWarp);
        
        
        
        % NEW: 
        % We need to obtain the information to populate this table from a json
        % file or somewhere. It will have defaults. Right now leave it hardcoded
        % for tests. 
        moriRois=["ATR_roi1_L", "ATR_roi2_L", ""; ...
                  "ATR_roi1_R", "ATR_roi2_R", ""; ...
                  "CST_roi1_L", "CST_roi2_L", ""; ...
                  "CST_roi1_R", "CST_roi2_R", ""; ...
                  "CGC_roi1_L", "CGC_roi2_L", ""; ...
                  "CGC_roi1_R", "CGC_roi2_R", ""; ...
                  "HCC_roi1_L", "HCC_roi2_L", ""; ...
                  "HCC_roi1_R", "HCC_roi2_R", ""; ...
                  "FP_R"      , "FP_L"      , ""; ...
                  "FA_L"      , "FA_R"      , ""; ...
                  "IFO_roi1_L", "IFO_roi2_L", ""; ...
                  "IFO_roi1_R", "IFO_roi2_R", ""; ...
                  "ILF_roi1_L", "ILF_roi2_L", ""; ...
                  "ILF_roi1_R", "ILF_roi2_R", ""; ...
                  "SLF_roi1_L", "SLF_roi2_L", ""; ...
                  "SLF_roi1_R", "SLF_roi2_R", ""; ...
                  "UNC_roi1_L", "UNC_roi2_L", ""; ...
                  "UNC_roi1_R", "UNC_roi2_R", ""; ...
                  "SLF_roi1_L", "SLFt_roi2_L", ""; ...
                  "SLF_roi1_R", "SLFt_roi2_R", ""; ...
                  "Left-LGN_dil-2", "V1_dil-1_L", ""; ...
                  "Right-LGN_dil-2", "V1_dil-2_R", ""] ;
        labels = ["Left Thalamic Radiation", ...
				  "Right Thalamic Radiation", ...
                  "Left Corticospinal", ...
                  "Right Corticospinal", ...
                  "Left Cingulum Cingulate", ...
                  "Right Cingulum Cingulate", ...
                  "Left Cingulum Hippocampus",...
                  "Right Cingulum Hippocampus", ...
                  "Callosum Forceps Major", ...
                  "Callosum Forceps Minor", ...
                  "Left IFOF",...
                  "Right IFOF", ...
                  "Left ILF",...
                  "Right ILF", ...
                  "Left SLF",...
                  "Right SLF", ...
                  "Left Uncinate",...
                  "Right Uncinate", ...
                  "Left Arcuate",...
                  "Right Arcuate", ...
                  "Left Optic Radiation V1", ...
                  "Right Optic Radiation V1"]';
        slabels = ["LTR","RTR", "LCST","RCST", ...
                  "LCC", "RCC", "LCH","RCH", ...
                  "CFMaj", "CFMin", "LIFOF","RIFO", ...
                  "LILF","RILF", "LSLF","RSLF", ...
                  "LUF","RUF", "LAF","RAF", ...
                  "LORV1", "RORV1"]'; 
        nhlabels = ["TR","TR", "CST","CST", ...
                    "CC", "CC", "CH","CH", ...
                    "CFMAJ", "CFMIN", "IFOF","IFO", ...
                    "ILF","ILF", "SLF","SLF", ...
                    "UF","UF", "AF","AF", ...
                    "ORV1", "ORV1"]'; 
        tracts           = table();
        tracts.roi1      = moriRois(:,1);
        tracts.extroi1   = [repmat(".nii.gz",[length(labels),1])];
        tracts.roi2      = moriRois(:,2);
        tracts.extroi2   = [repmat(".nii.gz",[length(labels),1])];
        tracts.roi3      = moriRois(:,3);
        tracts.extroi3   = [repmat(".nii.gz",[length(labels),1])];
		tracts.dilroi1   = [repmat("",[length(labels)-2,1]);"2";"2"];
		tracts.dilroi2   = [repmat("",[length(labels)-2,1]);"";""];
		tracts.dilroi3   = [repmat("",[length(labels)-2,1]);"";""];
        tracts.label     = labels;
        tracts.fgnum     = string([1:height(tracts)]');
        tracts.hemi      = repmat(["Left","Right"]',[length(labels)/2,1]);
        tracts.slabel    = slabels;
        tracts.shemi     = repmat(["L","R"]',[length(labels)/2,1]);
        tracts.nhlabel   = nhlabels;
        tracts.wbt       = [repmat(true,[length(labels)-4,1]); 0;0;0;0];
        tracts.usecortex = repmat(false,[length(labels),1]);
        tracts.quiet     = repmat("-quiet",[length(labels),1]);
        tracts.force     = repmat("-force",[length(labels),1]);
        tracts.maxlen    = repmat(200,[length(labels),1]);
        tracts.cmaxlen   = strcat("-maxlength ",num2str(tracts.maxlen));
        tracts.minlen    = repmat(10,[length(labels),1]);
        tracts.cminlen   = strcat("-minlength ",num2str(tracts.minlen));
        % Codify all options in the tract name? Or in a folder?
        tracts.fname     = strcat(tracts.slabel,".tck");
        tracts.cfname    = strcat(tracts.slabel,"_clean.tck");
        tracts.fdir      = repmat(mrtrixDir,[length(labels),1]);
        tracts.fpath     = strcat(tracts.fdir,filesep,tracts.fname);
        % Cleaning options
        tracts.cfpath    = strcat(tracts.fdir,filesep,tracts.cfname);
        tracts.clean     = repmat(true,[length(labels),1]); % [repmat(true,[length(labels)-2,1]);false;false];
        tracts.nfibers   = zeros(size(labels));
        tracts.cnfibers  = zeros(size(labels));
        tracts.maxDist   = 3*ones(size(labels));
        tracts.maxLen    = 3*ones(size(labels));
        tracts.numNodes  = 100*ones(size(labels));
        tracts.meanmed   = repmat("median",[length(labels),1]);
        tracts.maxIter   = [3*ones(size(labels,1)-2,1);1;1];
        % added my lmx
        tracts.angle     = repmat(45,[length(labels),1]);
        tracts.algorithm = repmat('iFoD2',[length(labels),1]); 
        tracts.cutoff    = repmat(0.1, [length(labels),1]);

        
        % Obtain the segmented tracts. 
        % In the new version:
        %   - Uses mrtrix tools for the tracking/segmentation
        %   - If requested, it will return cleaned fibers. 
        %   - It saves everything to disk. We can do a cleaning routine at the
        %     end if required. 
        %   - It will save clip to roi == 1 and 0 files to disk, but 
        
        [fg_classified, fg_clean, fg] = RTP_TractsGet(dtFile, afq, tracts);
        

    %% Compute Tract Profiles
    
    if true % AFQ_get(afq,'compute profiles',ii)
        fprintf('\nComputing Tract Profiles for subject %s',sub_dirs{ii});
        % Determine how much to weight each fiber's contribution to the
        % measurement at the tract core. Higher values mean steaper falloff
        fWeight = AFQ_get(afq,'fiber weighting');
        % By default Tract Profiles of diffusion properties will always be
        % calculated
        [fa,md,rd,ad,cl,vol,TractProfile]=AFQ_ComputeTractProperties(...
                                                fg_clean, ...
                                                dt, ...
                                                afq.params.numberOfNodes, ...
                                                afq.params.clip2rois, ...
                                                sub_dirs{ii}, ...
                                                fWeight, ...
                                                afq, ...
												tracts);
	
	
        % Parameterize the shape of each fiber group with calculations of
        % curvature and torsion at each point and add it to the tract profile
        [curv, tors, TractProfile] = AFQ_ParamaterizeTractShape(fg_classified, TractProfile);
        
        % Calculate the volume of each Tract Profile
        TractProfile = AFQ_TractProfileVolume(TractProfile);
        
        % Add values to the afq structure
        afq = AFQ_set(afq,'vals','subnum',ii,'fa',fa,'md',md,'rd',rd,...
            'ad',ad,'cl',cl,'curvature',curv,'torsion',tors,'volume',vol);
        
        % Add Tract Profiles to the afq structure
        afq = AFQ_set(afq,'tract profile','subnum',ii,TractProfile);
        
        % If any other images were supplied calculate a Tract Profile for that
        % parameter
        numimages = AFQ_get(afq, 'numimages');
        if numimages > 0
            for jj = 1:numimages
                % Read the image file
                image = niftiRead(afq.files.images(jj).path{ii});
                % Check image header
                if ~all(image.qto_xyz(:) == image.sto_xyz(:))
                   image = niftiCheckQto(image);
                end  
                % Resample image to match dwi resolution if desired
                if AFQ_get(afq,'imresample')
                    image = mrAnatResampleToNifti(image, fullfile(afq.sub_dirs{ii},'bin','b0.nii.gz'),[],[7 7 7 0 0 0]);
                end
                % Compute a Tract Profile for that image
                imagevals = AFQ_ComputeTractProperties(fg_classified, image, afq.params.numberOfNodes, afq.params.clip2rois, sub_dirs{ii}, fWeight, afq);
                % Add values to the afq structure
                afq = AFQ_set(afq,'vals','subnum',ii,afq.files.images(jj).name, imagevals);
                clear imagevals
            end
        end
    else
        fprintf('\nTract Profiles already computed for subject %s',sub_dirs{ii});
    end
    
    % Save each iteration of afq run if an output directory was defined
    if ~isempty(AFQ_get(afq,'outdir')) && exist(AFQ_get(afq,'outdir'),'dir')
        if ~isempty(AFQ_get(afq,'outname'))
            outname = fullfile(AFQ_get(afq,'outdir'),AFQ_get(afq,'outname'));
        else
            outname = fullfile(AFQ_get(afq,'outdir'),['afq_' date]);
        end
        save(outname,'afq');
    end
    
    % clear the files that were computed for this subject
    clear fg fg_classified TractProfile
end  % Ends runsubs

