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
multishell          = opts.track.multishell;
useACT              = opts.track.mrtrix_useACT;
faFodThresh         = opts.track.faFodThresh;
% Life
life_runLife        = opts.track.life_runLife  ; 
life_discretization = opts.track.life_discretization;
life_num_iterations = opts.track.life_num_iterations;
life_test           = opts.track.life_test;
life_saveOutput     = opts.track.life_saveOutput;
life_writePDB       = opts.track.life_writePDB;
% ET
ET_angleValues      = opts.track.ET_angleValues;
ET_maxlength        = opts.track.ET_maxlength;
ET_minlength        = opts.track.ET_minlength;
ET_numberFibers     = opts.track.ET_numberFibers;
ET_stepSizeMm       = opts.track.ET_stepSizeMm;

                                                       
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
% Edit GLU 03.2020: Now ET is always run, if we add just a pair of options, it will be equivalent of doing it only once. Big difference is that now we are combining angle and length in the tractotrams. 
% TODO: clean the helop above to reflect the new function.
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
if notDefined('ET_stepSizeMm'), ET_stepSizeMm = 999; end
if notDefined('ET_numberFibers'), ET_numberFibers = 400000; end
if notDefined('ET_angleValues'), ET_angleValues = [45, 25, 5]; end
if notDefined('ET_maxlength'), ET_maxlength = [100, 150, 200]; end
if notDefined('ET_minlength'), ET_minlength = 20; end


if mrtrixVersion ~= 3
    error('Mrtrix3 supported only')
end


% Variables to consider here
% --- algo: make it use the defaults per every case. No switch-case for this one
% --- multishell or not: it seems that at this point it doesn't matter, it
%     was important for the FoD calculation. Now we can use ACT or not. 

% This is yes/no options afterwards
% --- LiFE or not: do it or not over the previous output
% --- Save pdb or not: do it or not over the previous output

% Generate the optionals here
% They will be empty strings if nothing is going to be added and mrtrix
% results will be used. 
if faFodThresh == 999
    faFodThreshStr = '';
else
    faFodThreshStr = [' -cutoff ' num2str(faFodThresh) ' '];
end

if ET_minlength == 999
    ET_minlengthStr = '';
else
    ET_minlengthStr = [' -minlength ' num2str(ET_minlength) ' '];
end

if ET_stepSizeMm == 999
    ET_stepSizeMmStr = '';
 else
    ET_stepSizeMmStr = [' -step ' num2str(ET_stepSizeMm) ' '];
end

% Create the optional str that will be added to the tckgen call
optionalStr = [faFodThreshStr ET_minlengthStr ET_stepSizeMmStr];

% Generate the appropriate UNIX command string.
[~, pathstr] = strip_ext(files.csd);
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
        fgFileName{na}=['ET_fibs' num2str(ET_numberFibers) ...
                                '_angle' strrep(num2str(ET_angleValues(na)),'.','p') ...
                                '_maxlen' strrep(num2str(ET_maxlength(na)),'.','p') ...
                                '.tck'];
        fgFileNameWithDir{na}=fullfile(fileparts(tck_file), fgFileName{na});
        cmd_str = ['tckgen -force ' ...
                    '-algo ' algo ...
                    optionalStr ...
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

if life_runLife
    % "clean" the tractogram before using it further
    % First, obtain the dirNames the mainLife that we already had was using, so
    % that we do not change much in the first iteration

    mrtrixFolderParts  = split(files.csd, filesep);
    % Obtain the session name. This is usually the zip name if it has not
    % been edited. 
    mrtrixDir  = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-1)), filesep);
    dtiDir     = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-2)), filesep);
    sessionDir = dtiDir;
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
    fg.name = [NAME '_LiFE'];
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
