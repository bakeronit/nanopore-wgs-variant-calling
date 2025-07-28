#!/usr/bin/env bash
# properties = {properties}
# This script gets executed on a compute node on the cluster

echo -e "JOB ID\t$PBS_JOBID"
module load apptainer/1.4.1
{exec_job}
