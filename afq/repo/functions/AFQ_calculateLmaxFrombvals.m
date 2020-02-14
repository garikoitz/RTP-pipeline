function lmax = AFQ_calculateLmaxFrombvals(dt6)
    % calculate max lmax (6->2, 15->4, 28->6, etc..)
 
    % GLU 08.2018
    
    % Loading the dt file containing all the paths to the fiels we need.
    dt_info = load(dt6);
    [PATHSTR,NAME,EXT] = fileparts(dt_info.files.alignedDwBvals);
    dt6Parts = split(dt6, filesep);
    AnalysisDir = strjoin(dt6Parts(1:(length(dt6Parts)-2)), filesep);
    bvalues = dlmread(fullfile(AnalysisDir, [NAME EXT]));
    howMany = sum(bvalues > 100);
    
    lmax=0;
    while (lmax+3)*(lmax+4)/2 <= howMany
       lmax = lmax + 2;
    end
   

end