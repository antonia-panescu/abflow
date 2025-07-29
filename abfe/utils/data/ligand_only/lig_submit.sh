#!/bin/bash
#SBATCH --nodes=1
## The use of ntasks-per-socket instead of ntasks-per-node is *very* important
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=24:00:00
#SBATCH -J tak
#SBATCH -p g-ran
#SBATCH --array=1-33
#SBATCH --output=output_%x_%A_%a.out   # Output file name (includes array ID)

# Load the relevant GROMACS module
# Load the relevant GROMACS module
module add GROMACS/2023.1-foss-2022a-CUDA-11.7.0

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
			-p	../../topol_$leg.top \
			-c	../../lig_$leg.gro \
			-n	../../index.ndx \
			-r	../../lig_$leg.gro \
			-maxwarn 3 \
			-o	$leg.enmin.$number.tpr 


gmx mdrun -s -ntmpi 1  -pin on -pinstride 1 -deffnm $leg.enmin.$number


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
			-p	../../topol_$leg.top \
			-c	../enmin/$leg.enmin."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r	../enmin/$leg.enmin."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o	$leg.nvt.$number.tpr 

gmx mdrun -s -ntmpi 1  -pin on -pinstride 1 -deffnm $leg.nvt.$number 

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


########################## NPT_b         ######################

cd ../npt_b

gmx grompp  -f	$leg.npt_b.$number.mdp \
			-p	../../topol_$leg.top \
			-c	../nvt/$leg.nvt."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r	../nvt/$leg.nvt."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o	$leg.npt_b.$number.tpr  

gmx mdrun -s -ntmpi 1  -pin on -pinstride 1 -deffnm $leg.npt_b.$number 

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

############################# NPT_pr     #########################

cd ../npt_pr

gmx grompp  -f	$leg.npt_pr.$number.mdp \
			-p	../../topol_$leg.top \
			-c ../npt_b/$leg.npt_b."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r ../npt_b/$leg.npt_b."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o $leg.npt_pr.$number.tpr 

gmx mdrun -s -ntmpi 1  -pin on -pinstride 1 -deffnm $leg.npt_pr.$number 

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
			-p	../../topol_$leg.top \
			-c ../nvt/$leg.nvt."$number"_nojump_mol.gro \
			-n	../../index.ndx \
			-r ../nvt/$leg.nvt."$number"_nojump_mol.gro \
			-maxwarn 3 \
			-o $leg.prod.$number.tpr

gmx mdrun -s -ntmpi 1  -pin on -pinstride 1 -deffnm $leg.prod.$number -v -dhdl dhdl

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


