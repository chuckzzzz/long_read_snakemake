#!/bin/bash
#SBATCH -A b1042               # Allocation
#SBATCH -p genomicsguestA                # Queue
#SBATCH -t 48:00:00             # Walltime/duration of the job
#SBATCH --mem=128G
#SBATCH --cpus-per-task=14
#SBATCH --output=/projects/b1042/YueLab/zzhang/ont_test/job_log/mapping.%j.%N.txt
#SBATCH --error=/projects/b1042/YueLab/zzhang/ont_test/job_log/mapping.%j.%N.err
source /home/zzj4347/.bashrc
conda activate ont
cd /projects/b1042/YueLab/zzhang/ont_test
snakemake --cores all -p
