#!/bin/bash
#SBATCH --job-name=$JOBNAME
#SBATCH --nodes=$NO_OF_NODES
#SBATCH --output=out.%x.%A.%a.out
#SBATCH --tasks-per-node=$TASKS_PER_NODE
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --time=24:00:00
#SBATCH --account=e280-biggin 
#SBATCH --partition=standard
#SBATCH --qos=standard


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
# ENMIN
############################################################################################################


# Loop through each window and run the command
for i in "${!windows[@]}"; do
    window="${windows[$i]}"
    
    # Navigate to the appropriate directory
    cd "$window/enmin"
    
    # Define the necessary variables
    window_type="${window:0:-3}"
    window_number="${window: -2}"
    
    # Construct and run the command using mpiexec
    srun --cpu-bind=cores gmx_mpi grompp -f "${window_type}.enmin.${window_number}.mdp" \
           -c "../../complex_${window_type}.gro" \
           -r "../../complex_${window_type}.gro" \
           -n "../../index_${window_type}.ndx" \
           -p "../../complex_${window_type}.top" \
           -o "topol.tpr" \
           -maxwarn 2
    
    # Go back to the starting directory
    cd ../../
done

srun --cpu-bind=cores gmx_mpi mdrun -multidir $MULTIDIR_LIST_STR_ENMIN

############################################################################################################
# NVT
############################################################################################################


# Loop through each window and run the command
for i in "${!windows[@]}"; do
    window="${windows[$i]}"
    
    # Navigate to the appropriate directory
    cd "$window/nvt"
    
    # Define the necessary variables
    window_type="${window:0:-3}"
    window_number="${window: -2}"
    
    # Construct and run the command using mpiexec
    srun --cpu-bind=cores gmx_mpi grompp -f "${window_type}.nvt.${window_number}.mdp" \
           -c "../enmin/confout.gro" \
           -r "../enmin/confout.gro" \
           -n "../../index_${window_type}.ndx" \
           -p "../../complex_${window_type}.top" \
           -o "topol.tpr" \
           -maxwarn 2
    
    # Go back to the starting directory
    cd ../../
done

srun --cpu-bind=cores gmx_mpi mdrun -multidir $MULTIDIR_LIST_STR_NVT 

############################################################################################################
# NPT_B
############################################################################################################


# Loop through each window and run the command
for i in "${!windows[@]}"; do
    window="${windows[$i]}"
    
    # Navigate to the appropriate directory
    cd "$window/npt_b"
    
    # Define the necessary variables
    window_type="${window:0:-3}"
    window_number="${window: -2}"
    
    # Construct and run the command using mpiexec
    srun --cpu-bind=cores gmx_mpi grompp -f "${window_type}.npt_b.${window_number}.mdp" \
           -c "../nvt/confout.gro" \
           -r "../nvt/confout.gro" \
           -n "../../index_${window_type}.ndx" \
           -p "../../complex_${window_type}.top" \
           -o "topol.tpr" \
           -maxwarn 2
    
    # Go back to the starting directory
    cd ../../
done

srun --cpu-bind=cores gmx_mpi mdrun -multidir $MULTIDIR_LIST_STR_NPT_B 

############################################################################################################
# NPT_PR
############################################################################################################


# Loop through each window and run the command
for i in "${!windows[@]}"; do
    window="${windows[$i]}"
    
    # Navigate to the appropriate directory
    cd "$window/npt_pr"
    
    # Define the necessary variables
    window_type="${window:0:-3}"
    window_number="${window: -2}"
    
    # Construct and run the command using mpiexec
    srun --cpu-bind=cores gmx_mpi grompp -f "${window_type}.npt_pr.${window_number}.mdp" \
           -c "../npt_b/confout.gro" \
           -r "../npt_b/confout.gro" \
           -n "../../index_${window_type}.ndx" \
           -p "../../complex_${window_type}.top" \
           -o "topol.tpr" \
           -maxwarn 2
    
    # Go back to the starting directory
    cd ../../
done

srun --cpu-bind=cores gmx_mpi mdrun -multidir $MULTIDIR_LIST_STR_NPT_PR 

############################################################################################################
# PROD
############################################################################################################


# Loop through each window and run the command
for i in "${!windows[@]}"; do
    window="${windows[$i]}"
    
    # Navigate to the appropriate directory
    cd "$window/prod"
    
    # Define the necessary variables
    window_type="${window:0:-3}"
    window_number="${window: -2}"
    
    # Construct and run the command using mpiexec
    srun --cpu-bind=cores gmx_mpi grompp -f "${window_type}.prod.${window_number}.mdp" \
           -c "../npt_pr/confout.gro" \
           -r "../npt_pr/confout.gro" \
           -n "../../index_${window_type}.ndx" \
           -p "../../complex_${window_type}.top" \
           -o "topol.tpr" \
           -maxwarn 2
    
    # Go back to the starting directory
    cd ../../
done

srun --cpu-bind=cores gmx_mpi mdrun -multidir $MULTIDIR_LIST_STR_PROD -replex 200 -hrex -plumed ../../plumed.dat -dlb no
