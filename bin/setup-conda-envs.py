#!/usr/bin/env python

import os
from time import sleep

# addresses issue with parallel installation of conda envs: https://github.com/nextflow-io/nextflow/issues/1819

# meant to be run in genelab-utils conda env, inside the workflow root directory

## setting up a dictionary with yaml file and name with digest nextflow would create
    # these need to be updated if anything changes in the yaml files (like versions)
envs_dict = {}

envs_dict['trim-galore.yaml'] = "trim-galore-60fba8f2aafaec1bf1a22b398630243a"
envs_dict['QC.yaml'] = "QC-7853e05fc8ff88ec7805f4a89d7c42dd"
envs_dict['dp-tools.yaml'] = "dp-tools-35f84ab11c0f0fe15c56f66bcc2b15d6"
envs_dict['gtf-to-bed.yaml'] = "gtf-to-bed-182b867c918aee90de26b6d83492ef9a"
envs_dict['bismark.yaml'] = "bismark-82258a372f4dfaad72d87c5ef67fc5f5"
envs_dict['methylkit.yaml'] = "methylkit-33b71761ca76f60c9c031cda49844fd3"
envs_dict['align-qc.yaml'] = "align-qc-bdae5ff37e3903a19417a93877402765"
envs_dict['nugen-trim.yaml'] = "python2-a6e35c3da0172ba5bd255dfdbc24c0a4"

## relative path to conda environment yamls
env_yaml_path = "config/software/conda-envs/"

# checking this is being run from the expected location
if not os.path.isdir(env_yaml_path):

    print(f"\n  This script is expected to be run in the root directory of the workflow,")
    print(f"  but the {env_yaml_path} cannot be found.")
    print("\nExiting for now.\n")
    exit(1)

else:

    print("\nAttempting to setup conda environments...")

## getting environment variable for the conda prefix
CONDA_PREFIX = os.getenv('CONDA_PREFIX')

### running setup for each if they don't exist already
for curr_yaml in envs_dict.keys():

    curr_env = curr_yaml.split(".")[0]

    curr_location = os.path.join(CONDA_PREFIX, "envs", envs_dict[curr_yaml])

    curr_yaml_full_path = os.path.join(env_yaml_path, curr_yaml)

    # only trying to create if environment doesn't already exist
    if not os.path.isdir(curr_location):

        curr_command = f"mamba env create -p {curr_location} -f {curr_yaml_full_path}"

        print(f"\n\n  Setting up {curr_env} environment with:")
        print(f"    {curr_command}\n")
        sleep(1)

        os.system(curr_command)
    
    else:

        print(f"\n\n  {curr_env} already detected, run the following if wanting to remove it:")
        print(f"    conda env remove -p {curr_location}")
        sleep(1)

print("\n  All conda environments that will be used by the workflow are present.\n")
