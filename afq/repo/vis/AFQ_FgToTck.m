function AFQ_FgToTck(fg, pathName, fib2mt, createObj)
    for fiberNum  = 1:length(fg)
        myfg      = fg(fiberNum); 
        [fiberPath, fname, fext] = fileparts(pathName);
        % Add the name of .mat this tracts originated
        fiberName = [fname '_' strrep(myfg.name,' ','') '.tck'];
        objName   = [fname '_' strrep(myfg.name,' ','') '.obj'];
        % Usually tck files stored in mrtrix folder
        if fib2mt;fiberPath = strrep(fiberPath,'fibers','mrtrix'); end
        if ~exist(fiberPath,'dir'); mkdir(fiberPath);end
        AFQ_WriteMrtrixTck(myfg, fullfile(fiberPath, fiberName))
        if createObj;AFQ_fg2MNIObj(myfg,'fname', fullfile(fiberPath,objName));end
    end
end
