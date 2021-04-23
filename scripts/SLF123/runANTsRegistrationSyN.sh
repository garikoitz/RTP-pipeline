#!/bin/bash
module load $ants_ver

date;
printf "#### output dir: $output \n"
printf "#### template: $template \n"
printf "#### image: $img \n"
printf "#### ants version: $ants_ver \n"
printf "############################### \n"

cmd="antsRegistrationSyN.sh -d 3 -o $output/ants -f $template -m $img"

echo $cmd
eval $cmd
