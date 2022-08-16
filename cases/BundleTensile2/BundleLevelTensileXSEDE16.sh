#!/bin/bash

#SBATCH --partition=RM


# For shared memory parallelization
#SBATCH --cpus-per-task=16
#SBATCH	--ntasks=1
#SBATCH	--nodes=1

#SBATCH --time 0-01:00:00

#SBATCH	--account=watsonkh
#SBATCH	--mail-user=watsonkh@vcu.edu
#SBATCH	--mail-type=END
#SBATCH	--mail-type=BEGIN


#SBATCH --job-name=16Scale

#SBATCH --mem=100000

julia -t $SLURM_CPUS_PER_TASK ./cases/BundleTensile2/BundleLevelTensileXSEDE.jl ./cases/BundleTensile2/test_cases_5/high_align_high_branch/InputLongitudinalTensile.toml ./cases/BundleTensile2/test_cases_5/high_align_high_branch/