# SLF123 script package

This script package includes skull stripping, non-linear registration, transformation and probabilistic map generation for the roisList.txt for the subjectList.txt
database.

## Files

- config.json: set up all the variables, including the script you want to run:
  1. registration.sh: Non-linear registration to the MNI template (you can choose any other reference you want by editing the proper variables in config.json)
  2. skullStripping.sh: if subjects have skull, perform skull stripping before non-linear registration
  3. transform.sh: transform native ROIs to the warped brain to obtain ROIs in the reference space

- if ROIs are in mat, before transform  convert them to Nifti by running rois_mat2nifti.m
- probabilisticMap.m : once registration and transformation are completed, run this matlab script to obtain the probabilistic maps. Furtheremore, it contains the option of cleaning the maps, i.e., select one ROI for the voxels that are labeled by more than one ROI (it selects the ROI with higher number of subjects for those voxels).

- probabilisticMap folder: it contains all the matlab functions needed for the probabilisticMap.m

## Usage

```bash
python SLF123.py config.json
```
