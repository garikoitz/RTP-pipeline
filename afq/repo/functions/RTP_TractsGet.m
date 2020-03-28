% function [fg_classified,fg_unclassified,classification,fg]=RTP_TractsGet(...
function [fg_classified,fg_clean,fg]=RTP_TractsGet(dt6File,afq) 
% Categorizes each fiber in a group into one of the 20 tracts defined in
% the Mori white matter atlas. 
%
%  [fg_classified, fg_unclassified, classification, fg] = ...
%      AFQ_SegmentFiberGroups(dt6File, fg, [Atlas='MNI_JHU_tracts_prob.nii.gz'], ...
%      [useRoiBasedApproach=true], [useInterhemisphericSplit=true], [antsInvWarp]);
%
%  Fibers are segmented in two steps. Fibers become candidates for a fiber
%  group if the pass through the 2 waypoint ROIs that define the
%  tracjectory of the tract. Then each fiber is compared to a fiber
%  proability map and high probability fibers are retained in the group.
%  The segmentation alogrithm is based on:
%
%  Hua K, Zhang J, Wakana S, Jiang H, Li X, Reich DS, Calabresi PA, Pekar
%  JJ, van Zijl PC, Mori S. 2008. Tract probability maps in stereotaxic
%  spaces: analyses of white matter anatomy and tract-specific
%  quantification. Neuroimage 39(1):336-47.
%
%  Zhang, W., Olivi, A., Hertig, S. J., van Zijl, P., & Mori, S. (2008).
%  Automated fiber tracking of human brain white matter using diffusion
%  tensor imaging. Neuroimage, 42(2), 771-777.
%
% Input parameters:
% dt6File                  - Either the dt6 structure or a path to the
%                            dt6.mat file
% fg                       - A file with previously tracked elsewhere fibers
%                          to be categorized.
% Atlas                    - probabilistic atlas defining probabilities for
%                          each voxel to be passed by a fiber within each of
%                          atlas fiber groups. We usually use Mori atlas
%                          supplied with fsl: MNI_JHU_tracts_prob.nii.gz.
%                          This atlas is not symmetric. For a "symmetrified"
%                          atlas use 'MNI_JHU_tracts_prob_Symmetric.nii.gz'
%                          but we strongly recommend using the original
%                          atlas.
% useInterhemisphericSplit - cut fibers crossing between hemispheres with a
%                          midsaggital plane below z=-10. This is to get
%                          rid of CST fibers that appear to cross at the
%                          level of the brainstem
% useRoiBasedApproach      - use the approach describing in Zhang (2008) Neuroimage 42.
%                           For each of the 20 Mori Groups 2 critical ROIs
%                           are computed by spatially transforming ROIs
%                           provided in templates/MNI_JHU_tracts_ROIs. A
%                           fiber becomes a candidate to being labeled as a
%                           part of a given Mori group if this fiber
%                           "passes through" both critical ROIs for that
%                           Mori group. Our modification of Zhang (2008)
%                           approach: In case a single fiber is a candidate
%                           for >1 Mori group, respective cumulative
%                           probabilities are computed with probabilistic
%                           Mori atlas, then compared.
%                           useRoiBasedApproach can take the following
%                           values: (1) 'false' (to not use the approach);
%                           (2) 'true'(use the approach; the minimal
%                           distance from a fiber to ROI to count as "a
%                           fiber  is crossing the ROI" minDist=2mm); (3) a
%                           scalar value for minDist in mm; (4) a vector
%                           with the first value being minDist and the
%                           second value being a flag 1/0 for whether ROIs
%                           should be recomputed (and overwritten) for this
%                           subject (1, by default). This is useful because
%                           sometimes if you  are rerunning
%                           dtiFindMoriTracts with different parameters,
%                           you do not need to recompute ROIs.  E.g., to
%                           avoid recomputing  ROIs and use minDist of 4mm
%                           one would pass [useRoiBasedApproach=[4 0]];
%  antsInvWarp              - Spatial normalization computed with ANTS. If
%                           a path to a precomputed ANTS warp is passed in
%                           then it will be used to transform the MNI ROIs
%                           to native space
%
% Output parameters:
% fg_ classified  - fibers structure containing all fibers assigned to
%                   one of Mori Groups. Respective group labeles are stored
%                   in fg.subgroups field.
% fg_unclassified - fiber structure containing the rest of the (not Mori) fibers.
% classification  - This variable gives the fiber group names and the group
%                   fiber group number for each fiber in the input group
%                   fg.  classification is a structure with two fields.
%                   classification.names is a cell array where each cell is
%                   the name of that fiber group. For example
%                   classification.names{3} = 'Corticospinal tract L'.
%                   classification.index is a vector that defines which
%                   group number each fiber in the origional fiber group
%                   was assigned to. For example
%                   classification.index(150)=3 means that fg.fibers(150)
%                   is part of the corticospinal tract fiber group.  The
%                   values in classification may not match the origional
%                   fiber group because of pre-processing.  However they
%                   will match the output fg which is the origional group
%                   with preprocessing.
% fg              - This is the origional pre-segmented fiber group.  It
%                   may differ slightly from the input due to preprocessing
%                   (eg splitting fibers that cross at the pons, removing 
%                   fibers that are too short)
%                    
% Example:
%    AFQdata = '/home/jyeatman/matlab/svn/vistadata/AFQ';
%    dt6File = fullfile(AFQdata, 'subj2', 'dt6.mat');
%    fg      = AFQ_SegmentFiberGroups(dt6File);
%
% See also: dtiSplitInterhemisphericFibers
%
% (c) Vistalab

%% Check arguments
disp(['[RTP_TractsGet] Using MORI ROIs to select fibers from whole brain tractogram']);

%% Read the data
% Load the dt6 file
if ischar(dt6File)
    dt = dtiLoadDt6(dt6File);
    baseDir = fileparts(dt6File); 
else
    dt = dt6File;
    baseDir = fileparts(dt.dataFile);
    dt6File = dt.dataFile;
end


% Make sure fiber folder exist before tracking or checking
fibDir = fullfile(baseDir,'fibers');
if ~exist(fibDir,'dir'); mkdir(fibDir); end
% Make sure mrtrix folder exist before tracking or checking
mrtrixDir = fullfile(baseDir,'mrtrix');
if ~exist(mrtrixDir,'dir'); mkdir(mrtrixDir); end

% The Whole Brain Tactrography will only be done here. 
if sum(afq.tracts.wbt) > 0
	% Check if it exist and force is false:
	if exist(fullfile(fibDir,'WholeBrainFG.mat'),'file') && afq.force==false
		fprintf('WBT file %s exists and will not be overwritten\n',fullfile(fibDir,'WholeBrainFG.mat'))
        fg = AFQ_get(afq,'wholebrain fiber group',1);
	else
   		 fprintf('\nPerforming whole-brain tractograpy\n');
         %% Track with mrtrix if the right files are there
         mrtrixpaths   = AFQ_get(afq,'mrtrix paths',afq.currentsub);
         mrtrixVersion = 3;
         [status, results, fg, pathstr] = AFQ_mrtrix_track(mrtrixpaths, ...
                                                           mrtrixpaths.wm, ... % roi, (wm = wmMask)
                                                           mrtrixpaths.wm_dilated,...  % wmMask
                                                           afq.params.track.mrTrixAlgo, ...
                                                           afq.params.track.ET_numberFibers,...
                                                           [],...
                                                           [],...
                                                            1, ...
                                                            mrtrixVersion, ...
                                                            afq.params);
   		 % Save it
   		 dtiWriteFiberGroup(fg,fullfile(fibDir,'WholeBrainFG.mat'));
   		 % Set the path to the fibers in the afq structure
   		 afq = AFQ_set(afq,'wholebrain fg path','subnum',1,fullfile(fibDir,'WholeBrainFG.mat'));
	end
end

% Start simple, with the existing tracts that we only want to do the tckedit
afq.tracts.fname     = strcat(afq.tracts.slabel,".tck");
afq.tracts.cfname    = strcat(afq.tracts.slabel,"_clean.tck");
afq.tracts.fdir      = repmat(mrtrixDir,[height(afq.tracts),1]);
afq.tracts.fpath     = strcat(afq.tracts.fdir,filesep,afq.tracts.fname);
afq.tracts.cfpath    = strcat(afq.tracts.fdir,filesep,afq.tracts.cfname);
afq.tracts.nfibers   = zeros(height(afq.tracts),1);
afq.tracts.cnfibers  = zeros(height(afq.tracts),1);


tracts = afq.tracts;
for nt=1:height(tracts)

	% TODO: add the clipping, create fg_clip and write the fibers. 
    % Mount all the components of the call based on the options in tracts table
    ts     = tracts(nt,:);
    fprintf('[RTP_TractsGet] Getting %s ...\n', ts.label)
    if exist(ts.fpath,'file')
       fprintf('[RTP_TractsGet] Using exisging one because force=false\n')
    	% Read the tract, we want to have the same fg struct as before
        tract = fgRead(ts.fpath);
        if nt==1;fg_classified=tract;
        else; fg_classified(nt)=tract;end
	    if exist(ts.cfpath,'file') 
			clean_tract = fgRead(ts.cfpath);
	        % Add it to fg_clean
    	    if nt==1; fg_clean=clean_tract;
       		 else; fg_clean(nt)=clean_tract; end
        end
	else
	    % Add the path
        RoiPara = load(dt6File);
        fs_dir  = RoiPara.params.fs_dir;
        moridir = fullfile(fs_dir, 'ROIs');
        % Solve the dilate text
        if ts.dilroi1>0;dil1=strcat("_dil-",num2str(ts.dilroi1));else;dil1="";end 
        if ts.dilroi2>0;dil2=strcat("_dil-",num2str(ts.dilroi2));else;dil2="";end 
        if ts.dilroi3>0;dil3=strcat("_dil-",num2str(ts.dilroi3));else;dil3="";end 
        roi1    = fullfile(moridir, strcat(ts.roi1,dil1,ts.extroi1));
        roi2    = fullfile(moridir, strcat(ts.roi2,dil2,ts.extroi2));
        roi3    = "";
        if ~(strcmp(ts.roi3,"NO"))
            join("-include", fullfile(moridir, strcat(ts.roi3,dil3,ts.extroi3)));
        end
    
    
        % The most important is wbt (whole brain tractography), whether we want to
        % use it or track the tract ourselves (ex OR)
        if ts.wbt
        % TODO: This code is still working, make two wbt, one with interhemispheric fibers and one without 
        % Make an ROI for the mid saggital plane
	    % midSaggitalRoi= dtiRoiMakePlane([0, dt.bb(1, 2), dt.bb(1, 3); 0 , dt.bb(2, 2) , dt.bb(2, 3)], 'midsaggital', 'g');
	    % keep1=zeros(length(fg.fibers), size(moriRois, 1));
	    % keep2=zeros(length(fg.fibers), size(moriRois, 1));
	    % Find fibers that cross mid saggital plane
	    % [fgOut, contentiousFibers, InterHemisphericFibers] = dtiIntersectFibersWithRoi([], 'not', [], midSaggitalRoi, fg);
	    %NOTICE: ~keep3 (not "keep3") will mark fibers that DO NOT cross
	    %midSaggitalRoi.
	    % keep3=repmat(InterHemisphericFibers, [1 size(moriRois, 1)]);



       fprintf('\n[RTP_TractsGet] Using tckedit to separate %s\n', ts.label)
       % Select the tracts that go in to tckedit
       tracks_in = fullfile(char(ts.fdir),[fg.name '.tck']);
       cmd       = join(["tckedit", ...
                    "-include",roi1, "-include",roi2, roi3, ...
                    "-maxlength",ts.maxlen, "-minlength", ts.minlen, ...
                    tracks_in, ts.fpath]);
       spres     = AFQ_mrtrix_cmd(cmd);
       % Make it readable and writeable
       fileattrib(ts.fpath, '+w +x') 
       
    else
		% TODO: add better logic roi1, roi2, roi3, what is seed, what it is waypoint
        fprintf('\n[RTP_TractsGet] Using tckgen to create %s\n', ts.label)
        cmd       = join([ "tckgen", "-algorithm", ts.algorithm, ...
                            "-select 5000", ...
                            "-seed_image", roi1, ...
                            "-include", roi2, ...
                            "-seed_image", roi2, ...
                            "-include", roi1, ...
                            "-angle", ts.angle, "-cutoff", ts.cutoff, ...
                            "-minlength", ts.minlen, "-maxlength", ts.maxlen, ...
                            "-stop", afq.files.mrtrix.csd{1}, ts.fpath]);
        spres     = AFQ_mrtrix_cmd(cmd);
    end

    % Read the tract, we want to have the same fg struct as before
    tract = fgRead(ts.fpath);
    if nt==1; fg_classified=tract;
	else; fg_classified(nt)=tract; end

    % Update the value of the number of fibers
    ts.nfibers = size(tract.fibers,1);

    % If more than 10 fibers clean it, otherwise copy it as it is
    if ts.nfibers>10
       clean_tract = AFQ_removeFiberOutliers(tract,ts.maxDist,ts.maxLen,ts.numNodes,ts.meanmed,1,ts.maxIter);
	else
	   clean_tract = tract;
    end

	% Add it to fg_clean
    if nt==1; fg_clean=clean_tract;
    else; fg_clean(nt)=clean_tract; end

    AFQ_fgWrite(clean_tract, ts.cfpath,'tck');
    fileattrib(ts.cfpath, '+w +x') % make it readable and writeable
    % Update the value of the number of fibers
    ts.cnfibers = size(clean_tract.fibers,1);    
    
    % Update the table, maybe we updated some of the fields (e.g. nfibers)
    tracts(nt,:) = ts;
    fprintf('\n[RTP_TractsGet] ... done %s\n', ts.label)
end
end
afq = rmfield(afq,'tracts');
afq.tracts = tracts;


% Cut the fibers below the acpc plane, to disentangle CST & ATL crossing 
% at the level of the pons
% TODO: decide if we want to keep this
%{
if useInterhemisphericSplit
    fgname  = fg.name;
    fg      = dtiSplitInterhemisphericFibers(fg, dt, -10);
    fg.name = fgname;
end
%}
%{
  % Warp the fibers to the MNI standard space so they can be compared to the
  % template
  fg_sn = dtiXformFiberCoords(fg, invDef);
  
  % moriTracts.data is a an XxYxZx20 array contianing the 20 Mori probability
  % atlases (range is 0-100 where 100 represents p(1)).
  sz = size(moriTracts.data);
  % fg_sn fiber coords are in MNI space- now convert them to atlas space by
  % applying the affine xform from the atlas NIFTI header. Since the atlas is
  % already in MNI space, this transform will just account for any
  % translation and scale differences between the atlas maps and the MNI
  % template used to compute our sn.
  fgCoords = mrAnatXformCoords(moriTracts.qto_ijk, horzcat(fg_sn.fibers{:}));
  clear fg_sn;   % what we need from fg_sn is now stored in fgCoords
  fgLen = cellfun('size',fg.fibers,2);
  
  % Now loop over the 20 atlases and get the atlas probability score for each
  % fiber point. We collapse the scores across all points in a fiber by taking
  % the mean. Below, we will use these 20 mean scores to categorize the fibers.
  % TO DO: consider doing something more sophisticated than taking the mean.
  fp = zeros(sz(4),numel(fg.fibers));
  for(ii=1:sz(4))
      % Get the Mori atlas score for each point in the fibers using
      % trilinear interpolation.
      p = myCinterp3(double(moriTracts.data(:,:,:,ii))/100, sz([1,2]), sz(3), fgCoords(:,[2,1,3]));
      % The previous line interpolated one giant array with all fiber points
      % concatenated. The next loop will separate the coordinates back into
      % fibers and take the mean score for the points within each fiber.
      fiberCoord = 1;
      for(jj=1:numel(fg.fibers))
          fp(ii,jj) = nanmean(p([fiberCoord:fiberCoord+fgLen(jj)-1]));
          fiberCoord = fiberCoord+fgLen(jj);
      end
  end
  clear p fgCoords;
%}


%% Find fibers that pass through both waypoint ROIs eg. Wakana 2007
%Warp Mori ROIs to individual space; collect candidates for each fiber
%group based on protocol of 2 or > ROIs a fiber should travel through. The
%following ROIs are saved within
%trunk/mrDiffusion/templates/MNI_JHU_tracts_ROIs folder and are created
%using MNI template as described in Wakana et al.(2007) Neuroimage 36 with
%a single modification: For SLFt Roi2, they recommend drawing the ROI at
%the AC level, whereas we use a slice just inferior of CC splenium. The
%reason for this modification is that Wakana et al. ACPC aligned images
%appear different from MNI images (the latter we use for defininng ROIs).
%If defining SLFt-Roi2 on a slice actually at the AC level (althought
%highly consistently across human raters), many SLFt fibers were not
%correctly labeled as they extend laterally into temporal lobe just above
%the aforementioned ROI plane.


% Remove below too, we have done the tracking above. We can use the atlas above to clean the outlier fibers better than with the current method.

if false % useRoiBasedApproach 
    % A 20x2 cell array containing the names of both waypoint ROIs for each
    % of 20 fiber groups

    % Make an ROI for the mid saggital plane
    midSaggitalRoi= dtiRoiMakePlane([0, dt.bb(1, 2), dt.bb(1, 3); 0 , dt.bb(2, 2) , dt.bb(2, 3)], 'midsaggital', 'g');
    keep1=zeros(length(fg.fibers), size(moriRois, 1)); 
    keep2=zeros(length(fg.fibers), size(moriRois, 1));
    % Find fibers that cross mid saggital plane
    [fgOut, contentiousFibers, InterHemisphericFibers] = dtiIntersectFibersWithRoi([], 'not', [], midSaggitalRoi, fg); 
    %NOTICE: ~keep3 (not "keep3") will mark fibers that DO NOT cross
    %midSaggitalRoi.
    keep3=repmat(InterHemisphericFibers, [1 size(moriRois, 1)]);
    fgCopy=fg; fgCopy.subgroup=[];
    for nt=1:size(moriRois, 1)
        % Load the nifit image containing ROI-1 in MNI space
		% find ROI location
       RoiPara = load(dt6File);
       fs_dir = RoiPara.params.fs_dir;
       moridir = fullfile(fs_dir, 'MORI');
       ROI_img_file=fullfile(moridir, moriRois{nt, 1});
        % Transform ROI-1 to an individuals native space
%         if recomputeROIs
%             % Default is to use the spm normalization unless a superior
%             % ANTS normalization was passed in
%             if exist('antsInvWarp','var') && ~isempty(antsInvWarp)
%                 outfile = fullfile(fileparts(dt6File),'ROIs',moriRois{nt, 1});
%                 roi1(nt) = ANTS_CreateRoiFromMniNifti(ROI_img_file, antsInvWarp, [], outfile);
%             else
%                 [RoiFileName, invDef, roi1(nt)]=dtiCreateRoiFromMniNifti(dt6File, ROI_img_file, invDef, true);
%             end
%         else
%             RoiFileName=fullfile(fileparts(dt6File), 'ROIs',  [prefix(prefix(ROI_img_file, 'short'), 'short') '.mat']);
%            roi1(nt) = dtiReadRoi(RoiFileName);  
             RoiFileName=fullfile(fileparts(dt6File), 'ROIs',  [prefix(prefix(moriRois{nt, 1}, 'short'), 'short') '.mat']);
             roi1(nt) = dtiImportRoiFromNifti(ROI_img_file);
             if ~exist(fullfile(fileparts(dt6File), 'ROIs'), 'dir')
                 mkdir(fullfile(fileparts(dt6File), 'ROIs'))
             end
             dtiWriteRoi(roi1(nt), RoiFileName)
%        end
        % Find fibers that intersect the ROI
        [fgOut,contentiousFibers, keep1(:, nt)] = dtiIntersectFibersWithRoi([], 'and', minDist, roi1(nt), fg);
        keepID1=find(keep1(:, nt));
        % Load the nifit image containing ROI-2 in MNI space
         ROI_img_file=fullfile(moridir, ['MORI_', moriRois{nt, 2}]);

%         % Transform ROI-2 to an individuals native space
%         if recomputeROIs
%             [RoiFileName, invDef, roi]=dtiCreateRoiFromMniNifti(dt6File, ROI_img_file, invDef, true);
%             % Default is to use the spm normalization unless a superior
%             % ANTS normalization was passed in
%             if exist('antsInvWarp','var') && ~isempty(antsInvWarp)
%                 outfile = fullfile(fileparts(dt6File),'ROIs',moriRois{nt, 2});
%                 roi2(nt) = ANTS_CreateRoiFromMniNifti(ROI_img_file, antsInvWarp, [], outfile);
%             else
%                 [RoiFileName, invDef, roi2(nt)]=dtiCreateRoiFromMniNifti(dt6File, ROI_img_file, invDef, true);
%             end   
%         else
%             RoiFileName=fullfile(fileparts(dt6File), 'ROIs',  [prefix(prefix(ROI_img_file, 'short'), 'short') '.mat']);
%            roi2(nt) = dtiReadRoi(RoiFileName);
             roi2(nt) = dtiImportRoiFromNifti(ROI_img_file);
             RoiFileName=fullfile(fileparts(dt6File), 'ROIs',  [prefix(prefix(moriRois{nt, 2}, 'short'), 'short') '.mat']);
             dtiWriteRoi(roi1(nt), RoiFileName)
%        end
        %To speed up the function, we intersect with the second ROI not all the
        %fibers, but only those that passed first ROI.
        fgCopy.fibers=fg.fibers(keepID1(keepID1>0));
        [a,b, keep2given1] = dtiIntersectFibersWithRoi([], 'and', minDist, roi2(nt), fgCopy);
        keep2(keepID1(keep2given1), nt)=true; 
    end
    clear fgOut contentiousFibers keepID
    %Note: forceps major and minor should NOT have interhemipsheric fibers
    % excluded
    keep3(:, 9:10)=keep3(:, 9:10).*0;
    % fp is the variable containing each fibers score for matching the
    % atlas. We will set each fibers score to 0 if it does not pass through
    % the necessary ROIs
    fp=ones(size(keep3'));
   
    fp(~(keep1'&keep2'&~keep3'))=0;
    %Also note: Tracts that cross through slf_t rois should be automatically
    %classified as slf_t, without considering their probs.
    fp(19, (keep1(:, 19)'&keep2(:, 19)'&~keep3(:, 19)'))=max(fp(:));
    fp(20, (keep1(:, 20)'&keep2(:, 20)'&~keep3(:, 20)'))=max(fp(:));
    
end



% TODO: don't know who did this, remove it

%{
% temporarily
sz(4)=20;
%% Eliminate fibers that don't match any of the atlases very well
% We have a set of atlas scores for each each fiber. To categorize the
% fibers, we will find the atlas with the highest score (using 'sort').
[atlasScore,atlasInd] = sort(fp,1,'descend');
% Eliminate fibers that don't match any of the atlases very well:
unclassified=atlasScore(1,:)==0;
goodEnough = atlasScore(1,:)~=0;
curAtlasFibers = cell(1,sz(4));
for ii=1:sz(4)
    curAtlasFibers{ii} = find(atlasInd(1,:)==ii & goodEnough);
end

% We now have a cell array (curAtlasFibers) that contains 20 arrays, each
% listing the fiber indices for the corresponding atlas group. E.g.,
% curAtlasFibers{3} is a list of indices into fg.fibers that specify the
% fibers belonging to group 3.

%Create a FG for unclassified fibers
fg_unclassified=fg;
fg_unclassified.name=[fg.name ' not Mori Groups'];
fg_unclassified.fibers = fg.fibers(unclassified); %prepare fg output
if ~isempty(fg.seeds)
    fg_unclassified.seeds = fg.seeds(unclassified);
end
fg_unclassified.subgroup=zeros(size(fg.fibers(unclassified)))'+(1+sz(4));
fg_unclassified.subgroupNames(1)=struct('subgroupIndex', 1+sz(4), 'subgroupName', 'NotMori');

% Create a fiber group for classified fibers
fg_classified = fg;
% Modify fg.fibers to discard the fibers that didn't make it into any of
% the atlas groups:
fg_classified.name=[fg.name ' Mori Groups'];
fg_classified.fibers = fg.fibers([curAtlasFibers{:}]);
if ~isempty(fg.seeds)
    fg_classified.seeds = fg.seeds([curAtlasFibers{:}],:);
end
%}

% TODO: check if this is necessary or not
% fg_classified.subgroup = zeros(1,numel(fg_classified.fibers));
% fg_clean.subgroup = zeros(1,numel(fg_clean.fibers));

%{
% We changed the size of fg_classified.fibers by discarding the
% uncategorized fibers, so we need to create a new array to categorize the
% fibers. This time we make an array with one entry corresponding to each
% fiber, with integer values indicating to which atlas group the
% corresponding fiber belongs.
fg_classified.subgroup = zeros(1,numel(fg_classified.fibers));
curInd = 1;
for(ii=1:numel(curAtlasFibers))
    fg_classified.subgroup(curInd:curInd+numel(curAtlasFibers{ii})-1) = ii;
    %fghandle.subgroup(curInd:curInd+numel(curAtlasFibers{ii})-1) = ii;
    curInd = curInd+numel(curAtlasFibers{ii});
    %Save labels for the fiber subgroups within the file
    fg_classified.subgroupNames(ii)=struct('subgroupIndex', ii, 'subgroupName', labels(ii));
    %fghandle.subgroupNames(ii)=struct('subgroupIndex', ii, 'subgroupName', labels(ii));
end
% The number of fibers in the origional preprocessed fiber group (fg) is
% equalt to the number of fibers in fg_classified + fg_unclassified. To
% create a new structure denoting the classification of each fiber group in
% the origional structure first preallocate an array the length of the
% origional (preprocessed) fiber group.
fiberIndex = zeros(length(fg.fibers),1);

% Populate the structure denoting the fiber group number that each fiber in
% the original wholebrain group was assigned to
for ii = 1 : length(curAtlasFibers)
    fiberIndex(curAtlasFibers{ii}) = ii;
    names{ii} = fg_classified.subgroupNames(ii).subgroupName;
end

% Create a structure with fiber indices and group names.
classification.index = fiberIndex;
classification.names = names;

% Convert the fiber group into an array of fiber groups
fg_classified = fg2Array(fg_classified);

% Flip fibers so that each fiber in a fiber group passes through roi1
% before roi2
for ii = 1:size(moriRois, 1)
    fg_classified(ii) = AFQ_ReorientFibers(fg_classified(ii),roi1(ii),roi2(ii));
end
%}






return;






