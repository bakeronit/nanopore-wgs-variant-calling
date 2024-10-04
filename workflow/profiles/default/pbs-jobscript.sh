#!/usr/bin/env bash
# properties = {properties}
# This script gets executed on a compute node on the cluster

echo -e "JOB ID\t$PBS_JOBID"
echo "================================="

{exec_job}