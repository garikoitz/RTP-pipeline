#! /usr/bin/env python

# Parse a config file and create a dtiInit params json file.
def parse_config(input_file,
                 output_file,
                 input_dir,
                 output_dir,
                 bvec_dir,
                 bval_dir,
                 nifti_dir,
                 anat_dir,
                 fs_dir,
                 tractparams_dir):
    import os
    import json
    import glob
    import shutil

    if not os.path.isfile(input_file):
        manifest = "/flywheel/v0/manifest.json"
        input_file = manifest
        MANIFEST=True
    else:
        MANIFEST=False

    # Read the config json file
    with open(input_file, 'r') as jsonfile:
        config_json = json.load(jsonfile)

    if MANIFEST:
        print "Loading default configuration from %s" % input_file
        manifest_config = dict.fromkeys(config_json['config'].keys())
        for k in manifest_config.iterkeys():
            manifest_config[k] = config_json['config'][k]['default']
        config = dict()
        config['params'] = manifest_config
    else:
        # Rename the config key to params
        print "Parsing %s" % input_file
        config = dict()
        if config_json['config'].has_key('config'):
            config['params'] = config_json['config']['config']
        else:
            config['params'] = config_json['config']

    # Handle the 'track' fields
    config['params']['track'] = {}
    config['params']['track']['faMaskThresh']        = config['params']['track_faMaskThresh']
    config['params']['track']['faFodThresh']         = config['params']['track_faFodThresh']

    config['params']['track']['mrtrix_useACT']       = config['params']['mrtrix_useACT']
    config['params']['track']['mrtrix_autolmax']     = config['params']['mrtrix_autolmax']
    config['params']['track']['mrtrix_lmax']         = config['params']['mrtrix_lmax']
    config['params']['track']['mrTrixAlgo']          = config['params']['mrtrix_mrTrixAlgo']

    config['params']['track']['get_vofparc']         = config['params']['get_vofparc']

    config['params']['track']['sift_runSift']        = config['params']['sift_runSift']
    config['params']['track']['sift_nFibers']        = config['params']['sift_nFibers']

    config['params']['track']['life_runLife']        = config['params']['life_runLife']
    config['params']['track']['life_discretization'] = config['params']['life_discretization']
    config['params']['track']['life_num_iterations'] = config['params']['life_num_iterations']
    config['params']['track']['life_test']           = config['params']['life_test']
    config['params']['track']['life_saveOutput']     = config['params']['life_saveOutput']
    config['params']['track']['life_writePDB']       = config['params']['life_writePDB']

    config['params']['track']['ET_numberFibers']     = config['params']['ET_numberFibers']
    if MANIFEST:
        config['params']['track']['ET_angleValues']      = [ float(x) for x in config['params']['ET_angleValues']['default'].split(',') ]
        config['params']['track']['ET_maxlength']        = [ float(x) for x in config['params']['ET_maxlength']['default'].split(',') ]
    else:
        config['params']['track']['ET_angleValues']      = [ float(x) for x in config['params']['ET_angleValues'].split(',') ]
        config['params']['track']['ET_maxlength']        = [ float(x) for x in config['params']['ET_maxlength'].split(',') ]

    config['params']['track']['ET_minlength']        = config['params']['ET_minlength']
    config['params']['track']['ET_stepSizeMm']       = config['params']['ET_track_stepSizeMm']

    # Remove the other track_ fields
    del config['params']['track_faMaskThresh']
    del config['params']['track_faFodThresh']

    del config['params']['mrtrix_useACT']
    del config['params']['mrtrix_autolmax']
    del config['params']['mrtrix_lmax']
    del config['params']['mrtrix_mrTrixAlgo']

    del config['params']['get_vofparc']

    del config['params']['sift_runSift']
    del config['params']['sift_nFibers']

    del config['params']['life_runLife']
    del config['params']['life_discretization']
    del config['params']['life_num_iterations']
    del config['params']['life_test']
    del config['params']['life_saveOutput']
    del config['params']['life_writePDB']

    del config['params']['ET_numberFibers']
    del config['params']['ET_angleValues']
    del config['params']['ET_maxlength']
    del config['params']['ET_minlength']
    del config['params']['ET_track_stepSizeMm']

    # Add input directory for dtiInit
    config['input_dir']       = input_dir
    config['output_dir']      = output_dir
    config['bvec_dir']        = bvec_dir
    config['bval_dir']        = bval_dir
    config['nifti_dir']       = nifti_dir
    config['anat_dir']        = anat_dir
    config['fs_dir']          = fs_dir
    config['tractparams_dir'] = tractparams_dir
    # new parameters: save_output (delcare if save .zip) and segmentSLF
    # (declare if apply segmentSLF)
    config['params']['save_output']     = config['params']['save_output']
   # del config['params']['save_output']
    config['params']['segmentSLF']      = config['params']['segmentSLF']
   # del config['params']['segmentSLF']
    # Add additional keys
    config['params']['run_mode'] = [],
    config['params']['outdir'] = []
    config['params']['input_dir'] = input_dir
    config['params']['output_dir'] = output_dir

    # Write out the modified configuration
    with open(output_file, 'w') as config_json:
        json.dump(config, config_json, sort_keys=True, indent=4, separators=(',', ': '))

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--input_file', default='/flwywheel/v0/config.json', help='Full path to the input file.')
    ap.add_argument('--output_file', default='/flywheel/v0/json/params.json', help='Full path to the output file.')
    ap.add_argument('--input_dir', default='/flwywheel/v0/input', help='Full path to the input file.')
    ap.add_argument('--output_dir', default='/flywheel/v0/output', help='Full path to the output file.')
    ap.add_argument('--bvec_dir', default='/flwywheel/v0/input/bvec', help='Full path to the input file.')
    ap.add_argument('--bval_dir', default='/flywheel/v0/input/bval', help='Full path to the output file.')
    ap.add_argument('--nifti_dir', default='/flwywheel/v0/input/dwi', help='Full path to the input file.')
    ap.add_argument('--anat_dir', default='/flywheel/v0/input/anat', help='Full path to the output file.')
    ap.add_argument('--fs_dir', default='/flywheel/v0/input/fs', help='Full path to the freesurfer directory.')
    ap.add_argument('--tractparams_dir', default='/flywheel/v0/input/tractparams', help='Full path to the parameters file.')

    args = ap.parse_args()

    parse_config(args.input_file,
                 args.output_file,
                 args.input_dir,
                 args.output_dir,
                 args.bvec_dir,
                 args.bval_dir,
                 args.nifti_dir,
                 args.anat_dir,
                 args.fs_dir,
                 args.tractparams_dir)


    # Example, for copy and pase
   # python parse_config.py --input_file /data/localhome/glerma/soft/RTP-pipeline/example_config.json \
   #         --output_file /data/localhome/glerma/soft/RTP-pipeline/example_output_parsed_kkkkkkk.json \
   #         --input_dir  /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/input \
   #         --output_dir  /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/output \
   #         --bvec_dir  /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/input/bvec \
   #         --bval_dir  /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/input/bval \
   #         --nifti_dir  /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/input/dwi \
   #         --anat_dir  /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/input/anat \
   #         --fs_dir  /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/input/fs \
   #         --tractparams_dir /black/localhome/glerma/TESTDATA/FS/17_CAMINO_6835_docker/pipeline/input/tractparams
