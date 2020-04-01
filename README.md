[![Docker Pulls](https://img.shields.io/docker/pulls/scitran/afq-pipeline.svg)](https://hub.docker.com/r/scitran/afq-pipeline/)
[![Docker Stars](https://img.shields.io/docker/stars/scitran/afq-pipeline.svg)](https://hub.docker.com/r/scitran/afq-pipeline/)

# RTP pipeline
Reproducible Tract Profiles pipeline for diffusion imaging. 

Reproducible Tract Profiles (RTP) comprises a set of methods to manage and analyze diffusion weighted imaging (DWI) data for reproducible tractography. The tools take MRI data from the scanner and process them through a series of analysis implemented as Docker containers that are integrated into a modern neuroinformatics platform (Flywheel, Docker, Singularity). The platform guarantees that the entire pipeline can be re-executed, using the same data and computational parameters. We describe the DWI analysis tools that are used to identify the positions of a user-defined number of tracts and their diffusion profiles. The combination of these three components defines a system that transforms raw data into reproducible tract profiles for publication.

Please [check the wiki](https://github.com/vistalab/RTP-pipeline/wiki) for the most updated installation and usage manuals. 
