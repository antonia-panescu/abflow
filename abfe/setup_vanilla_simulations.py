import os
import shutil
import re
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.coordinates.GRO import GROWriter

from abfe.utils.utils import merge_gro_files, add_water2topology, addtoppar2top, copy_jobscripts_vanilla
from abfe.utils.extract_membrane import extract_membrane
from abfe.utils.smirnoff_param_func import smirnoff_parameterize
from abfe.utils.gromacs_helpers_single_chain import *


class SimulationSetup:
    def __init__(self, base_path, charmm_folder, protein_folder, sdf_folder,
                 crystal_water_gro, mdp_templates_path, jobscript_template_path,
                 solvate_water_count=10286, crystal_water_count=136, num_nodes=2):
        self.base_path = Path(base_path)
        self.charmm_folder = Path(charmm_folder)
        self.protein_folder = Path(protein_folder)
        self.sdf_folder = Path(sdf_folder)
        self.crystal_water_gro = Path(crystal_water_gro)
        self.mdp_templates_path = Path(mdp_templates_path)
        self.jobscript_template_path = Path(jobscript_template_path)
        self.solvate_water_count = solvate_water_count
        self.crystal_water_count = crystal_water_count
        self.num_nodes = num_nodes

    def setup_simulation(self, compound_name, replicate: int = 1):
        complex_folder = self.base_path / f"A2A_{compound_name}"
        
        # name each repeat folder “vanilla_rep_N”
        vanilla_name = f"vanilla_rep_{replicate}"
        vanilla_path = complex_folder / vanilla_name
        vanilla_path.mkdir(parents=True, exist_ok=True)

        sdf_file = self.sdf_folder / f"{compound_name}.sdf"
        shutil.copy(sdf_file, vanilla_path)

        dest_charmm = vanilla_path / self.charmm_folder.name
        shutil.copytree(self.charmm_folder, dest_charmm)

        print(f"Setup complete in {vanilla_path}")
        self.run_vanilla_setup(vanilla_path, compound_name)

    def run_vanilla_setup(self, vanilla_path, compound_name):
        os.chdir(vanilla_path)
        sdf_file = f"{compound_name}.sdf"

        # Helper functions
        param_lig(sdf_file)
        self.copy_parameterized_protein()
        organize_topologies()
        addligtop(sdf_file)
        addliggro(sdf_file)
        addligposres(sdf_file)
        extract_membrane_from_charmm()
        copy_complex2membrane()
        addmemtop(sdf_file)
        solvate_system(self.solvate_water_count)  
        cp_mdps(self.mdp_templates_path)
        merge_gro_files('solv_fix.gro', self.crystal_water_gro)
        add_water2topology(self.crystal_water_count)
        addions('solv_fix_crystal_water.gro')
        create_index('system_solv_ions.gro')
        addtoppar2top("topol.top")
        copy_jobscripts_vanilla(self.jobscript_template_path, compound_name, self.num_nodes)
        os.chdir('..')

    def copy_parameterized_protein(self):
        for item in self.protein_folder.iterdir():
            if item.is_file():
                shutil.copy(item, item.name)
            elif item.is_dir():
                shutil.copytree(item, item.name)

    def setup_all(self, repeats: int = 1):
        """
        Loop over every ligand in `sdf_folder`, running `repeats` independent setups.

        Args:
            repeats: how many vanilla_rep_N folders per ligand to produce.
        """
        compound_list = [p.stem for p in self.sdf_folder.glob("*.sdf")]
        for compound in compound_list:
            for rep in range(1, repeats + 1):
                print(f"Setting up {compound}, replica {rep}")
                self.setup_simulation(compound, replicate=rep)
