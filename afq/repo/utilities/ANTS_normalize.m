function ANTS_normalize(image, template, outname)
% Normalize an image to the MNI template using ANTS

% Use the new MNI template by default
if notDefined('template') 
    tdir = fullfile(AFQ_directories,'templates','mni_icbm152_nlin_asym_09a_nifti');
    template = fullfile(tdir,'mni_icbm152_t2_tal_nlin_asym_09a.nii');
end
outname = prefix(prefix(image))
% Normalize the image to the template with ANTS
cmd = ['antsIntroduction.sh '...
       '-d 3 '...
       '-r ' template ' '...
       '-i ' image ' '...
       '-o ' outname];
system(cmd)


% GLU notes on the commands
% -d:  ImageDimension: 2 or 3 (for 2 or 3 Dimensional registration)
% -r:  Reference image
% -i:  Input image
% -o:  OUTPREFIX; A prefix that is prepended to all output files.





% Use Nonlinear antsRegistrationSyN.sh for registration of diffusion with
% anatomical
image    = 'dwi.nii.gz';
anat     = 't1w_acpc.nii.gz';
outname  = 'dwi2t1';
cmd = ['antsRegistrationSyNQuick.sh '...
       '-d 3 '...
       '-f ' anat ' '...
       '-m ' image ' '...
       '-o ' outname ' '...
       '-n 6'];
system(cmd)