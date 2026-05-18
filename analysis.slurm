#!/bin/bash
#SBATCH --job-name=dv_7
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cluster=gpu
#SBATCH --partition=eschneider
#SBATCH --mail-user=hjl28@pitt.edu 
#SBATCH --mail-type=END,FAIL
#SBATCH --time=8:00:00
#SBATCH --output=output.txt
#SBATCH --output=slurm_outputs/slurm-%j.out       
#SBATCH --error=slurm_outputs/slurm-%j.err      

module purge
module load power9_env/gcc-10.1.0_openmpi-4.0.5_cuda-11.1.0
module load hdf5
module load python/anaconda3-2019.10

python pre_plotting.py