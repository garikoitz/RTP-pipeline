#!/bin/bash
# Builds the gear/container
# The container can be exported using the export.sh script
# cd afq/source/bin
# . ./compile.sh
# cd ../../../dtiinit/source/
# . ./compile.sh
# cd ../../
# git add .
# git commit -m "Commiting before building $GEAR:$1"
GEAR=garikoitz/rtp-pipeline
docker build --tag $GEAR:$1 .
