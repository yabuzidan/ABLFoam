#!/bin/bash

#SBATCH --job-name="!@#"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1900M

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

############################################

# Set environment variables
module purge
module load mpi/openmpi
source /home/ubuntu/OpenFOAM/openfoam8/etc/bashrc

# Run script
time ./Allrun-parallel