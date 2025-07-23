# ABFE_Package
Code for setting up and running replica ABFE simulations

## HPC-specific dependencies

In addition to the Python environment (`abfe_environment.yml`), this workflow requires the following modules to be loaded on the HPC environment before running:

```bash
module load gromacs/2020.3-AVX2-GPU
module load plumed/2.8.1_cuda_mpi
```

## Usage

`SimulationSetup` class for preparing vanilla systems.

Example:
```python
from abfe.setup_vanilla_simulations import SimulationSetup

setup = SimulationSetup(
    base_path           = "/ABS/PATH/TO/WORK/FOLDER",                 # e.g. "/home/user/projects/a2a_run"
    charmm_folder       = "/ABS/PATH/TO/CHARMM_GMX",                  # …/system_setup/charmm/charmm_gmx
    protein_folder      = "/ABS/PATH/TO/PARAMETERISED_PROTEIN",       # …/protein_prep/protein_param
    sdf_folder          = "/ABS/PATH/TO/LIGAND_SDFS",                 # folder full of *.sdf
    crystal_water_gro   = "/ABS/PATH/TO/crystal_waters.gro",          # the file you merge in
    mdp_templates_path  = "/ABS/PATH/TO/MDP_TEMPLATES",               # dir that contains your .mdp files
    jobscript_template_path = "/ABS/PATH/TO/JOBSCRIPT_TEMPLATES",     # dir with SLURM/Archer scripts

    # the three below are optional (defaults shown):
    solvate_water_count = 10286,
    crystal_water_count = 136,
    num_nodes           = 2
)

# prepare every ligand in the SDF folder
setup.setup_all()

# or, for a single ligand, single run:
# setup.setup_simulation("cmp_123")

# all ligands, 5 repeats
# setup.setup_all(repeats=5)

# single ligand, custom repeats (three replicas for only cmp_42
# for replica in range(1, 4):
#    setup.setup_simulation("cmp_42", replicate=replica)
```

`ABFESetup` class for setting-up ABFE folders for simulation.

Example:
```from abfe.abfe_setup import ABFESetup

setup = ABFESetup(
    base_path="/path/to/simulations",
    ligands=["Ligand1", "Ligand2"],
    num_replicates=3,
    template_script_path="/path/to/job_template.sh",
    contd_script_path="/path/to/job_template_contd.sh",
    archer_nodes=22
)

setup.setup()
```



`ABFEAnalyzer` class for analyzing ABFE productions.

Example:
```from abfe_package.abfe_analyzer import ABFEAnalyzer

analyzer = ABFEAnalyzer(
    project_root="/path/to/simulations",
    abfe_subdirs=["abfe_van1_hrex_r1", "abfe_van2_hrex_r1", "abfe_van3_hrex_r1"],
)

# Run analysis:
# - Complex analysis for all folders/replicates
# - Ligand analysis only for the first replicate
# - Symlink ligand results for other replicates
# - Compile final results for all replicates

analyzer.run()
```
