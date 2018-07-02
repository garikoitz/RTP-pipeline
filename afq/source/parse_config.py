#! /usr/bin/env python

# Parse a config file and create a dtiInit params json file.
def parse_config(input_file, output_file, input_dir, output_dir):
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
    config['params']['track']['algorithm']          = config['params']['track_algorithm']
    config['params']['track']['angleThresh']        = config['params']['track_angleThresh']
    config['params']['track']['faMaskThresh']       = config['params']['track_faMaskThresh']
    config['params']['track']['faThresh']           = config['params']['track_faThresh']
    config['params']['track']['lengthThreshMm']     = [config['params']['track_minLengthThreshMm'], config['params']['track_maxLengthThreshMm']]
    config['params']['track']['nfibers']            = config['params']['track_nfibers']
    config['params']['track']['offsetJitter']       = config['params']['track_offsetJitter']
    config['params']['track']['seedVoxelOffsets']   = [config['params']['track_seedVoxelOffset_1'], config['params']['track_seedVoxelOffset_2']]
    config['params']['track']['stepSizeMm']         = config['params']['track_stepSizeMm']
    config['params']['track']['wPuncture']          = config['params']['track_wPuncture']
    config['params']['track']['whichAlgorithm']     = config['params']['track_whichAlgorithm']
    config['params']['track']['whichInterp']        = config['params']['track_whichInterp']

    # # TODO: Add and pull this from the manifest
    config['params']['track']['mrTrixAlgo']          = config['params']['mrtrix_mrTrixAlgo']
    config['params']['track']['multishell']          = config['params']['mrtrix_multishell']
    config['params']['track']['tool']                = config['params']['mrtrix_tool']
    config['params']['track']['life_runLife']        = config['params']['life_runLife']
    config['params']['track']['life_discretization'] = config['params']['life_discretization']
    config['params']['track']['life_num_iterations'] = config['params']['life_num_iterations']
    config['params']['track']['life_test']           = config['params']['life_test']
    config['params']['track']['life_saveOuput']      = config['params']['life_saveOuput']


    # Remove the other track_ fields
    del config['params']['track_algorithm']
    del config['params']['track_angleThresh']
    del config['params']['track_faMaskThresh']
    del config['params']['track_faThresh']
    del config['params']['track_maxLengthThreshMm']
    del config['params']['track_minLengthThreshMm']
    del config['params']['track_nfibers']
    del config['params']['track_offsetJitter']
    del config['params']['track_seedVoxelOffset_1']
    del config['params']['track_seedVoxelOffset_2']
    del config['params']['track_stepSizeMm']
    del config['params']['track_wPuncture']
    del config['params']['track_whichAlgorithm']
    del config['params']['track_whichInterp']
    del config['params']['mrtrix_mrTrixAlgo']
    del config['params']['mrtrix_multishell']
    del config['params']['mrtrix_tool']
    del config['params']['life_runLife']
    del config['params']['life_discretization']
    del config['params']['life_num_iterations']
    del config['params']['life_test']
    del config['params']['life_saveOuput']

    # Handle cutoffLower and cutoffUpper
    config['params']['cutoff'] = [config['params']['cutoffLower'], config['params']['cutoffUpper'] ]

    # Remove cutoff fields
    del config['params']['cutoffUpper']
    del config['params']['cutoffLower']

    # Add input directory for dtiInit
    config['input_dir'] = input_dir
    config['output_dir'] = output_dir

    # Add additional keys
    config['params']['run_mode'] = [],
    config['params']['outdir'] = []
    config['params']['outname'] = config['params']['AFQ_Output_Name']
    del config['params']['AFQ_Output_Name']
    config['params']['input_dir'] = input_dir
    config['params']['output_dir'] = output_dir

    # Write out the modified configuration
    with open(output_file, 'w') as config_json:
        json.dump(config, config_json, sort_keys=True, indent=4, separators=(',', ': '))

if __name__ == '__main__':

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--input_file', default='/flwywheel/v0/config.json', help='Full path to the input file.')
    ap.add_argument('--output_file', default='/flywheel/v0/json/dtiinit_params.json', help='Full path to the output file.')
    ap.add_argument('--input_dir', default='/flwywheel/v0/input', help='Full path to the input file.')
    ap.add_argument('--output_dir', default='/flywheel/v0/output', help='Full path to the output file.')
    args = ap.parse_args()

    parse_config(args.input_file, args.output_file, args.input_dir, args.output_dir)
