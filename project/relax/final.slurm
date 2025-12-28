#!/bin/bash
#SBATCH --job-name=relax
#SBATCH --nodes=10
#SBATCH --ntasks=32
#SBATCH --array=0-4
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --error=relax.err

module load python/3.7.0
module load spack/0.16.3
. /home/apps/spack/share/spack/setup-env.sh
spack load /vns4j3z    # your QE 7.0 hash



SYSTEMS=(
  g_scf
  zno_slab_scf
  mb_scf
  g_zno_scf
  mb_gzno
)


NAME=${SYSTEMS[$SLURM_ARRAY_TASK_ID]}


mkdir -p outputs scratch/$NAME logs

srun pw.x -in input/$NAME.in > outputs/$NAME.pw.out

