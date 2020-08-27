function [afq afq_C2ROI] = AFQ_run(sub_dirs, sub_group, afq)
% function [afq patient_data control_data norms abn abnTracts afq_C2ROI] = AFQ_run(sub_dirs, sub_group, afq)
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
    afq    = AFQ_set(afq,'current subject',ii);
    dtFile = fullfile(sub_dirs{ii},'dt6.mat');
    [fg_classified, fg_clean, fg, fg_C2ROI] = RTP_TractsGet(dtFile, afq);
    % Create new afq for the C2ROI tracts
    afq_C2ROI = afq;




   %% Compute Tracts
   % Load Dt6 File
   dt     = dtiLoadDt6(dtFile);
   fprintf('\n[AFQ_run] Computing Tract Profiles for subject %s',sub_dirs{ii});
   % Determine how much to weight each fiber's contribution to the
   % measurement at the tract core. Higher values mean stepper falloff
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
                                            afq);
   
   
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



    %% Now do the same for the afq_C2ROI
   [fa,md,rd,ad,cl,vol,TractProfile]=AFQ_ComputeTractProperties(...
                                           fg_C2ROI, ...
                                           dt, ...
                                           afq_C2ROI.params.numberOfNodes, ...
                                           afq_C2ROI.params.clip2rois, ...
                                           sub_dirs{ii}, ...
                                           fWeight, ...
                                           afq_C2ROI);
   
   
   % Parameterize the shape of each fiber group with calculations of
   % curvature and torsion at each point and add it to the tract profile
   [curv, tors, TractProfile] = AFQ_ParamaterizeTractShape(fg_C2ROI, TractProfile);
   
   % Calculate the volume of each Tract Profile
   TractProfile = AFQ_TractProfileVolume(TractProfile);
   
   % Add values to the afq structure
   afq_C2ROI = AFQ_set(afq_C2ROI,'vals','subnum',ii,'fa',fa,'md',md,'rd',rd,...
       'ad',ad,'cl',cl,'curvature',curv,'torsion',tors,'volume',vol);
   
   % Add Tract Profiles to the afq structure
   afq_C2ROI = AFQ_set(afq_C2ROI,'tract profile','subnum',ii,TractProfile);
 






















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
    
   %  % Save each iteration of afq run if an output directory was defined
   %  if ~isempty(AFQ_get(afq,'outdir')) && exist(AFQ_get(afq,'outdir'),'dir')
   %      if ~isempty(AFQ_get(afq,'outname'))
   %          outname = fullfile(AFQ_get(afq,'outdir'),AFQ_get(afq,'outname'));
   %      else
   %          outname = fullfile(AFQ_get(afq,'outdir'),['afq_' date]);
   %      end
   %      save(outname,'afq');
   %  end
    
    % clear the files that were computed for this subject
    clear fg fg_classified TractProfile
end  % Ends runsubs

