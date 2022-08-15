#!/bin/bash

#SBATCH --partition=basic


# For shared memory parallelization
#SBATCH --cpus-per-task=200
#SBATCH	--ntasks=1
#SBATCH	--nodes=1

#SBATCH --time 99-00:00:00

#SBATCH	--account=watsonkh
#SBATCH	--mail-user=watsonkh@vcu.edu
#SBATCH	--mail-type=END
#SBATCH	--mail-type=BEGIN


#SBATCH --job-name=HAHB-MergedSplits

#SBATCH --mem=100000


julia ./cases/BundleTensile2/BundleLevelTensileTransverse.jl ./cases/BundleTensile2/test_cases_5/high_align_high_branch/InputTransverseTensile.toml ./cases/BundleTensile2/test_cases_5/high_align_high_branch/