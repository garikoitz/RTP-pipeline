function [status, results, fg, pathstr] = AFQ_mrtrix_track( files, ...
                                                            roi, ...
                                                            mask, ...
                                                            algo, ...
                                                            nSeeds, ...
                                                            bkgrnd, ...
                                                            verbose, ...
                                                            clobber, ...
                                                            mrtrixVersion, ...
                                                            opts)
 % Fix this                                                           

multishell          = opts.multishell;
useACT              = opts.mrtrix_useACT;
faThresh            = opts.faThresh;
% Life
life_runLife        = opts.life_runLife  ; 
life_discretization = opts.life_discretization;
life_num_iterations = opts.life_num_iterations;
life_test           = opts.life_test;
life_saveOutput     = opts.life_saveOutput;
life_writePDB       = opts.life_writePDB;
% ET
ET_runET            = opts.ET_runET;
ET_numberFibers     = opts.ET_numberFibers;
ET_angleValues      = opts.ET_angleValues;
ET_minlength        = opts.ET_minlength;
ET_maxlength        = opts.ET_maxlength;

                                                        
                                                       
%
% function [status, results, fg, pathstr] = mrtrix_track(files, roi, mask, algo, nSeeds, bkgrnd, verbose)
%
% Provided a csd estimate, generate estimates of the fibers starting in roi
% and terminating when they reach the boundary of mask
%
% Parameters
% ----------
% files: structure, containing the filenames generated using mrtrix_init.m
% roi: string, filename for a .mif format file containing the ROI in which
%      to place the seeds. Use the *_wm.mif file for Whole-Brain
%      tractography.
% mask: string, filename for a .mif format of a mask. Use the *_wm.mif file for Whole-Brain
%      tractography.
% algo: Tracking algorithm: it was 'mode' before. Specify it directly in afq.param.track.mrTrixAlgo
% nSeeds: The number of fibers to generate.
% bkgrnd: on unix, whether to perform the operation in another process
% verbose: whether to display standard output to the command window.
% clobber: Whether or not to overwrite the fiber group if it was already
%          computed
%
% Franco, Bob & Ariel (c) Vistalab Stanford University 2013
% Edit GLU 06.2016 added mrTrix versioning
% Mode now has been modified with algo, and it is set in the afq structure
% directly in AFQ_Create. Be careful, the options for mrTrix2 and mrTrix3 are
% very different. See AFQ_Create for all the available options.
% Edit GLU 06.2018 Added Life: after we get the tck we run it through the LiFE pipeline, so that 
% the fiber cleaning is done in origin. Then it is the cleaned WholeBrain
% converted to pdb and continues the normal afq pipeline.
% Edit GLU 08.2018: added Ensemble Tractography

status = 0; results = [];
if notDefined('verbose'),  verbose = false;end
if notDefined('bkgrnd'),    bkgrnd = false;end
if notDefined('clobber'),  clobber = false;end
if notDefined('mrtrixVersion'),    mrtrixVersion = 3; end
if notDefined('multishell'),  multishell = false;end
if notDefined('useACT'),  useACT = false;end

% LiFE
if notDefined('life_runLife'), life_runLife = true; end
if notDefined('life_discretization'), life_discretization = 360; end
if notDefined('life_num_iterations'), life_num_iterations = 100; end
if notDefined('life_test'), life_test = false; end
if notDefined('life_saveOutput'), life_saveOutput = false; end
if notDefined('life_writePDB'), life_writePDB = false; end
% Ensemble Tractography
if notDefined('ET_runET'), ET_runET = true; end
if notDefined('ET_numberFibers'), ET_numberFibers = 200000; end
if notDefined('ET_angleValues'), ET_angleValues = [47.2, 23.1, 11.5, 5.7, 2.9]; end
if notDefined('ET_minlength'), ET_minlength = 20; end
if notDefined('ET_maxlength'), ET_maxlength = 250; end


if mrtrixVersion ~= 3
    error('Mrtrix3 supported only')
end


% Variables to consider here
% --- algo: make it use the defaults per every case. No switch-case for this one
% --- multishell or not: it seems that at this point it doesn't matter, it
%     was important for the FoD calculation. Now we can use ACT or not. 

% This is 2x2 options
% --- ET or not: this generates 2 options per each case
% --- ACT or not: this generates 2 options per each case

% This is yes/no options afterwards
% --- LiFE or not: do it or not over the previous output
% --- Save pdb or not: do it or not over the previous output

% Generate the optionals here
% They will be empty strings if nothing is going to be added and mrtrix
% results will be used. 
if faThresh == 999
    faThreshStr = '';
else
    faThreshStr = [' -cutoff ' num2str(faThresh) ' '];
end

if ET_minlength == 999
    ET_minlengthStr = '';
else
    ET_minlengthStr = [' -minlength ' num2str(ET_minlength) ' '];
end

if ET_maxlength == 999
    ET_maxlengthStr = '';
else
    %ET_maxlengthStr = [' -maxlength ' num2str(ET_maxlength) ' '];
end

%optionalStr = [faThreshStr ET_minlengthStr ET_maxlengthStr];
optionalStr = [faThreshStr ET_minlengthStr];

% Generate the appropriate UNIX command string.
[~, pathstr] = strip_ext(files.csd);
if ET_runET
    if useACT
        disp('Running Ensemble Tractography with mrTrix3 and ACT.');
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_', algo, ...
                                          '-',num2str(nSeeds),'_ET_ACT.tck'));
        numconcatenate = [];
        for na=1:length(ET_angleValues)
            fgFileName{na}=['ET_fibs' num2str(ET_numberFibers) '_angle' strrep(num2str(ET_angleValues(na)),'.','p') '_ACT.tck'];
            fgFileNameWithDir{na}=fullfile(fileparts(tck_file), fgFileName{na});
            cmd_str = ['tckgen -force ' ...
                          '-algo ' algo optionalStr ' ' ...
                          '-backtrack -crop_at_gmwmi -info ' ...
                          '-seed_gmwmi ' files.gmwmi ' ' ...
                          '-act ' files.tt5 ' ' ...
                          '-angle ' num2str(ET_angleValues(na)) ' ' ...
                          '-select ' num2str(ET_numberFibers) ' ' ...
                          '-maxlength ' num2str(ET_maxlength(na)) ' ' ...
                          files.csd ' ' ...
                          fgFileNameWithDir{na}];
            % Run it, if the file is not there (this is for debugging)
            if ~exist(fgFileNameWithDir{na},'file')
                [status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);
            end
            numconcatenate = [numconcatenate, ET_numberFibers];
        end
        fg = et_concatenateconnectomes(fgFileNameWithDir, tck_file, numconcatenate, 'tck'); 
    else
        disp('Running Ensemble Tractography with mrTrix3 and no ACT.');
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_', algo, ...
                                          '-',num2str(nSeeds),'_ET.tck'));
        numconcatenate = [];
        for na=1:length(ET_angleValues)
            fgFileName{na}=['ET_fibs' num2str(ET_numberFibers) '_angle' strrep(num2str(ET_angleValues(na)),'.','p') '.tck'];
            fgFileNameWithDir{na}=fullfile(fileparts(tck_file), fgFileName{na});
            cmd_str = ['tckgen -force ' ...
                        '-algo ' algo optionalStr ' ' ...
                        '-seed_image ' roi ' ' ...
                        '-mask ' mask ' ' ...
                        '-angle ' num2str(ET_angleValues(na)) ' ' ...
                        '-select ' num2str(ET_numberFibers) ' ' ...
                        '-maxlength ' num2str(ET_maxlength(na)) ' ' ...
                        files.csd ' ' ...
                        fgFileNameWithDir{na}];
            % Run it, if the file is not there (this is for debugging)
            if ~exist(fgFileNameWithDir{na},'file')
                [status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);
            end
            numconcatenate = [numconcatenate, ET_numberFibers];
        end
        fg = et_concatenateconnectomes(fgFileNameWithDir, tck_file, numconcatenate, 'tck'); 

    end
else
    if useACT
        disp('Running tractography with no ET, with mrTrix3 and with ACT.');
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_', algo, ...
                                          '-',num2str(nSeeds),'_ACT.tck'));
        cmd_str = ['tckgen -force ' ...
                          '-algo ' algo ' ' ...
                          '-backtrack -crop_at_gmwmi -info ' ...
                          '-seed_gmwmi ' files.gmwmi ' ' ...
                          '-act ' files.tt5 ' ' ...
                          '-select ' num2str(nSeeds) ' ' ...
                          files.csd ' ' ...
                          tck_file];
    else
        disp('Running tractography with no ET, with mrTrix3 and with no ACT.');
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_', algo, ...
                                          '-',num2str(nSeeds),'.tck'));
        cmd_str = ['tckgen -force ' ...
                          '-algo ' algo optionalStr ' ' ...
                          '-seed_image ' roi ' ' ...
                          '-mask ' mask ' ' ...
                          '-select ' num2str(nSeeds) ' ' ...
                          files.csd ' ' ...
                          tck_file];
    end
    % Launch the resulting command 
    if ~(exist(tck_file,'file') ==2)  || clobber == 1
        [status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose, mrtrixVersion);
    else
        fprintf('\nFound fiber tract file: %s.\n Loading it rather than retracking', tck_file)
    end
    % Read the tck to the fg variable. 
    % If it is not changed below, it will be passed to next steps
    fg = fgRead(tck_file);
end

if life_runLife
    % "clean" the tractogram before using it further
    % First, obtain the dirNames the mainLife that we already had was using, so
    % that we do not change much in the first iteration

    mrtrixFolderParts  = split(files.csd, filesep);
    % Obtain the session name. This is usually the zip name if it has not
    % been edited. 
    mrtrixDir  = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-1)), filesep);
    dtiDir     = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-2)), filesep);
    sessionDir = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-3)), filesep);
    lifedir    = fullfile(dtiDir, 'LiFE');

    config.dtiinit             = dtiDir;
    config.track               = tck_file;
    config.life_discretization = life_discretization;
    config.num_iterations      = life_num_iterations;
    config.test                = life_test;
    % Change dir to LIFEDIR so that it writes everything there
    if ~exist(lifedir); mkdir(lifedir); end;
    cd(lifedir)

    disp('loading dt6.mat')
    disp(['Looking for file: ' fullfile(config.dtiinit, 'dt6.mat')])
    dt6 = load(fullfile(config.dtiinit, 'dt6.mat'))
    [~,NAME,EXT] = fileparts(dt6.files.alignedDwRaw);
    aligned_dwi = fullfile(sessionDir, [NAME,EXT])

    [ fe, out ] = life(config, aligned_dwi);

    out.stats.input_tracks = length(fe.fg.fibers);
    out.stats.non0_tracks = length(find(fe.life.fit.weights > 0));
    fprintf('number of original tracks	: %d\n', out.stats.input_tracks);
    fprintf('number of non-0 weight tracks	: %d (%f)\n', ...
             out.stats.non0_tracks, out.stats.non0_tracks / out.stats.input_tracks*100);
 
    if life_saveOutput
        disp('writing outputs')
        save('LiFE_fe.mat' ,'fe' , '-v7.3');
        save('LiFE_out.mat','out', '-v7.3');
    else
        disp('User selected not to write LiFE output')
    end

    % This is what we want to pass around
    fg = out.life.fg;
    % And I think I would need to write and substitute the non cleaned ET
    % tractogram tck with the new one...
    % Write file
    [PATHSTR,NAME,EXT] = fileparts(tck_file);
    tck_file =fullfile(PATHSTR, [NAME '_LiFE' EXT]);
    AFQ_fgWrite(fg, tck_file, 'tck');
end

if life_writePDB
    % This is the final output. Decide if we need the pdb output. 
    % It was removed because it was requiring huge amounts of RAM.
    % It can break an otherwise working gear. 
    % In any case, this should not affect the output, we want to pass fg to parent
    % function
    % The variable will still be called life_writePDB, though...
    % Convert the .tck fibers created by mrtrix to mrDiffusion/Quench format (pdb):
    % We will write both, but we want the cleaned ones to be used by vOF or any
    % other downstream code
    pdb_file = fullfile(pathstr,strcat(strip_ext(tck_file), '.pdb'));
    mrtrix_tck2pdb(tck_file, pdb_file);
end

end
