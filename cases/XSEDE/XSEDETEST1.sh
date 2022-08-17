#!/bin/bash

#SBATCH --partition=RM


# For shared memory parallelization
#SBATCH --cpus-per-task=1
#SBATCH	--ntasks=1
#SBATCH	--nodes=1

#SBATCH --time 0-01:00:00

#SBATCH	--account=mch210012p

#SBATCH --job-name=1ScaleTest

#SBATCH --mem=100000

julia -t $SLURM_CPUS_PER_TASK ./cases/XSEDE/XSEDETEST.jl ./cases/XSEDE/XSEDETEST.toml
