function [convert2RAS,dsline] = rtp_convert2RAScheck(fname)
%RTP_CONVERT2RASCHECK Summary of this function goes here
%   Detailed explanation goes here

% fname = '/Users/glerma/soft/rtp-pipeline/local/DTIINIT/input/T1w.nii.gz';

    convert2RAS     = false;
    [status,result] = AFQ_mrtrix_cmd(sprintf('mrinfo %s',fname));
    if status ~= 0 
        error('Mrtrixs mrinfo did return expected results, check if it installed properly')
    end
    
    TextAsCells     = regexp(result, '\n', 'split');
    dsline          = TextAsCells(contains(TextAsCells, "Data strides:"));
    dsline          = dsline{1};
    negs            = regexp(dsline, '-','once');
    if ~isempty(negs)
        disp('Negative direction detected, requires converting to RAS')
        convert2RAS=true;
    end
    dsline          = strrep(dsline,'Data strides:','');
    dsline          = strrep(dsline,'[','');
    dsline          = strrep(dsline,']','');
    dsline          = strrep(dsline,' ','');
    if (strcmp(dsline,'123') || strcmp(dsline,'1234'));
        disp('This file is RAS already, no conversion required')
    else    
        disp('This file is not RAS, conversion required')
        convert2RAS=true;
    end

end

