#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 01:24:52 2021

@author: llecca
"""

import argparse
import os
import json
import subprocess as sp

parser = argparse.ArgumentParser(description='''ANTsRegistrationSyN.py 'pathTo/config.json' ''')
# # Required positional argument
parser.add_argument('configFile', type=str, help='path to the config file')

args = parser.parse_args()
print('Read config file: ')
print(args.configFile)

# read the variables from json file:
with open(args.configFile,'r') as v:
    vars=json.load(v)

basedir=vars["config"]["basedir"]
codedir =vars["config"]["codedir"]
MNIpackage=vars["config"]["MNIpackage"]
MNItemplate=vars["config"]["MNItemplate"]
host =vars["config"]["host"]
qsub=vars["config"]["qsub"]
logdir=vars["config"]["logdir"]
script=vars["config"]["script"]

if not os.path.isdir(logdir): os.mkdir(logdir)

if host == "BCBL":
    mem=vars["BCBL"]["mem"]
    que=vars["BCBL"]["que"]
    core=vars["BCBL"]["core"]
    ants_ver=vars["BCBL"]["ants_ver"]
    fsl_ver=vars["BCBL"]["fsl_ver"]
    gcc_ver=vars["BCBL"]["gcc_ver"]
    mrtrix_ver=vars["BCBL"]["mrtrix_ver"]
elif host == "DIPC":
    mem=vars["DIPC"]["mem"]
    que=vars["DIPC"]["que"]
    core=vars["DIPC"]["core"]
    ants_ver=vars["DIPC"]["ants_ver"]
    fsl_ver=vars["DIPC"]["fsl_ver"]
    gcc_ver=vars["DIPC"]["gcc_ver"]
    mrtrix_ver=vars["DIPC"]["mrtrix_ver"]
    
os.chdir(codedir)

# run now the registration script
if script == "registration.sh":
    cmdstr = (f"{codedir}/registration.sh "+
               f"-b {basedir} " +
               f"-d {codedir} " +
               f"-p {MNIpackage} " +
               f"-t {MNItemplate} " +
               f"-h {host} " +
               f"-u {qsub} " +
               f"-g {logdir} " +
               f"-m {mem} " +
               f"-q {que} " + 
               f"-c {core} " +
               f"-a {ants_ver} ")
    print(cmdstr)
    sp.call(cmdstr, shell=True)
    
elif script == "transform.sh":
    cmdstr = (f"{codedir}/transform.sh "+
               f"-b {basedir} " +
               f"-d {codedir} " +
               f"-p {MNIpackage} " +
               f"-t {MNItemplate} " +
               f"-h {host} " +
               f"-u {qsub} " +
               f"-g {logdir} " +
               f"-m {mem} " +
               f"-q {que} " + 
               f"-c {core} " +
               f"-a {ants_ver} ")
    print(cmdstr)
    sp.call(cmdstr, shell=True)
    
elif script == "skullStripping.sh":
    cmdstr = (f"{codedir}/skullStripping.sh "+
               f"-b {basedir} " +
               f"-d {codedir} " +
               f"-p {MNIpackage} " +
               f"-t {MNItemplate} " +
               f"-h {host} " +
               f"-u {qsub} " +
               f"-g {logdir} " +
               f"-m {mem} " +
               f"-q {que} " + 
               f"-c {core} " +
               f"-a {ants_ver} " +
               f"-r {mrtrix_ver} " +
               f"-x {gcc_ver} " +
               f"-f {fsl_ver} ")
    print(cmdstr)
    sp.call(cmdstr, shell=True)

