#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --nodes=NO_OF_NODES
#SBATCH --output=out.%x.%A.%a.out
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=e280-biggin 
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --array=1-32

# Load the relevant GROMACS module
# Load the relevant GROMACS module
module load gromacs/2022.4

# to be careful we set the number of OpenMP threads before calling gromacs
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
string=$(sed -n "$SLURM_ARRAY_TASK_ID{p;q}" simulations_list.txt)

echo $string

leg="$(cut -d '.' -f1 <<< $string)"
number="$(cut -d'.' -f2 <<< $string)"


######################### ENMIN ###################################

cd ./$leg.$number
cd ./enmin

gmx grompp  -f	$leg.enmin.$number.mdp \
			-p	../../lig_$leg.top \
			-c	../../lig_$leg.gro \
			-n	../../index.ndx \
			-r	../../lig_$leg.gro \
			-maxwarn 3 \
			-o	$leg.enmin.$number.tpr 


srun --cpu-bind=cores gmx_mpi mdrun -s -ntmpi 1   -v -deffnm $leg.enmin.$number


gmx trjconv -f $leg.enmin.$number.gro \
			-s $leg.enmin.$number.tpr \
			-pbc nojump \
			-o $leg.enmin."$number"_nojump.gro \
			<< HERE
0
HERE
gmx trjconv -f $leg.enmin."$number"_nojump.gro \
			-s $leg.enmin.$number.tpr \
			-center -pbc mol \
			-o $leg.enmin."$number"_nojump_mol.gro \
			<< HERE
1
0
HERE

########################## Equilibration ######################
############################# NVT           ############################

cd ../nvt

gmx grompp  -f	$leg.nvt.$number.mdp \
			-p	../../lig_$leg.top \
			-c	../enmin/$leg.enmin."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r	../enmin/$leg.enmin."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o	$leg.nvt.$number.tpr 

srun --cpu-bind=cores gmx_mpi mdrun -s -ntmpi 1   -v -deffnm $leg.nvt.$number 

gmx trjconv -f $leg.nvt.$number.gro \
			-s $leg.nvt.$number.tpr \
			-pbc nojump \
			-o $leg.nvt."$number"_nojump.gro \
			<< HERE
0
HERE
gmx trjconv -f $leg.nvt."$number"_nojump.gro \
			-s $leg.nvt.$number.tpr \
			-center -pbc mol \
			-o $leg.nvt."$number"_nojump_mol.gro \
			<< HERE
1
0
HERE


########################## npt_B         ######################

cd ../npt_b

gmx grompp  -f	$leg.npt_b.$number.mdp \
			-p	../../lig_$leg.top \
			-c	../nvt/$leg.nvt."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r	../nvt/$leg.nvt."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o	$leg.npt_b.$number.tpr  

srun --cpu-bind=cores gmx_mpi mdrun -s -ntmpi 1   -v -deffnm $leg.npt_b.$number 

gmx trjconv -f $leg.npt_b.$number.gro \
			-s $leg.npt_b.$number.tpr \
			-pbc nojump \
			-o $leg.npt_b."$number"_nojump.gro \
			<< HERE
0
HERE
gmx trjconv -f $leg.npt_b."$number"_nojump.gro \
			-s $leg.npt_b.$number.tpr \
			-center -pbc mol \
			-o $leg.npt_b."$number"_nojump_mol.gro \
			<< HERE
1
0
HERE

############################# npt_PR     #########################

cd ../npt_pr

gmx grompp  -f	$leg.npt_pr.$number.mdp \
			-p	../../lig_$leg.top \
			-c ../npt_b/$leg.npt_b."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r ../npt_b/$leg.npt_b."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o $leg.npt_pr.$number.tpr 

srun --cpu-bind=cores gmx_mpi mdrun -s -ntmpi 1   -v -deffnm $leg.npt_pr.$number 

gmx trjconv -f $leg.npt_pr.$number.gro \
			-s $leg.npt_pr.$number.tpr \
			-pbc nojump \
			-o $leg.npt_pr."$number"_nojump.gro \
			<< HERE
0
HERE
gmx trjconv -f $leg.npt_pr."$number"_nojump.gro \
			-s $leg.npt_pr.$number.tpr \
			-center -pbc mol \
			-o $leg.npt_pr."$number"_nojump_mol.gro \
			<< HERE
1
0
HERE


############################# PROD        ############################

cd ../prod

gmx grompp  -f	$leg.prod.$number.mdp \
			-p	../../lig_$leg.top \
			-c ../npt_pr/$leg.npt_pr."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r ../npt_pr/$leg.npt_pr."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o $leg.prod.$number.tpr

srun --cpu-bind=cores gmx_mpi mdrun -s -ntmpi 1   -v -deffnm $leg.prod.$number -dhdl dhdl

gmx trjconv -f $leg.prod.$number.gro \
			-s $leg.prod.$number.tpr \
			-pbc nojump \
			-o $leg.prod."$number"_nojump.gro \
			<< HERE
0
HERE
gmx trjconv -f $leg.prod."$number"_nojump.gro \
			-s $leg.prod.$number.tpr \
			-center -pbc mol \
			-o $leg.prod."$number"_nojump_mol.gro \
			<< HERE
1
0
HERE
