#!/bin/bash
#SBATCH --job-name=$JOBNAME
#SBATCH --nodes=$NO_OF_NODES
#SBATCH --output=out.%x.%A.%a.out
#SBATCH --tasks-per-node=$TASKS_PER_NODE
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --time=24:00:00
#SBATCH --account=e280-biggin 
#SBATCH --partition=standard
#SBATCH --qos=$QOS


# to be careful we set the number of OpenMP threads before calling gromacs
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
# Disable backups
export GMX_MAXBACKUP=-1
export BASE_DIR=/work/e280/e280-Biggin/reub0138/gromacs_plumed_installation/
export GROMACS_DIR=${BASE_DIR}/gromacs-2020.6/
export PLUMED_DIR=${BASE_DIR}/plumed2-2.8.2/
export PATH=$PATH:${PLUMED_DIR}/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PLUMED_DIR}/lib
export PATH=$PATH:${GROMACS_DIR}/2020.6+plumed/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GROMACSD_DIR}/2020.6+plumed/lib


# Define the windows array
windows=($WINDOW_LIST_STR)



############################################################################################################
# PROD
############################################################################################################

srun --cpu-bind=cores gmx_mpi mdrun -multidir $MULTIDIR_LIST_STR_PROD -replex 200 -hrex -plumed ../../plumed.dat -dlb no -cpi state.cpt
