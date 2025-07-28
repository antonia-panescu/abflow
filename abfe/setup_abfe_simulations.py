import os
import logging
import sys
from pathlib import Path
from typing import List
from importlib.resources import files

from abfe.utils.abfe_helpers_hrex import (
    generate_boresch_restraints,
    copy_complex,
    create_fep_system,
    create_index,
    create_directories_with_MDP_pertuball_smaller_dumps,
    create_simulations_list,
    gen_hrex_submission_script
)


class ABFESetup:
    """
    ABFESetup handles the preparation of folders for ABFE simulations.
    """

    def __init__(self, base_path: str, ligands: List[str], num_replicates: int, archer_nodes: int = 22):
        self.base_path = Path(base_path).resolve()
        self.ligands = ligands
        self.num_replicates = num_replicates
        self.archer_nodes = archer_nodes

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('setup_abfe.log'),
                logging.StreamHandler()
            ]
        )
        logging.info(f"Initialized ABFESetup with base_path: {self.base_path}")

    def setup(self):
        for ligand in self.ligands:
            ligand_dir = self.base_path / ligand
            if not ligand_dir.exists():
                logging.warning(f"Ligand directory not found: {ligand_dir}")
                continue

            logging.info(f"Setting up ABFE folders for ligand: {ligand}")

            for rep in range(1, self.num_replicates + 1):
                suffix = f"van1_hrex_r{rep}"
                vanilla_folder_name = f"vanilla_rep_{rep}"
                self._setup_abfe_folder(ligand_dir, suffix, vanilla_folder_name)

    def _setup_abfe_folder(self, ligand_dir: Path, suffix: str, vanilla_folder_name: str):
        try:
            ligname = ligand_dir.name
            abfe_folder = ligand_dir / f"abfe_{suffix}"

            if abfe_folder.exists():
                logging.warning(f"{abfe_folder} already exists. Skipping...")
                return

            abfe_folder.mkdir(exist_ok=True)
            complex_path = abfe_folder / "complex"
            complex_path.mkdir(exist_ok=True)
            logging.info(f"Created folder: {abfe_folder} and {complex_path}")
            os.chdir(complex_path)  # << all operations now run inside complex/


            vanilla_path = ligand_dir / vanilla_folder_name
            if not vanilla_path.exists():
                logging.error(f"Vanilla directory does not exist: {vanilla_path}")
                return

            generate_boresch_restraints(vanilla_folder_name=vanilla_path, abfe_path=abfe_folder)

            copy_complex(vanilla_folder_name=vanilla_path, abfe_path=abfe_folder)
            create_fep_system(ligname, vanilla_path)
            create_index("complex_coul.gro", "index.ndx")
            protein_membrane_ligand_index = "Protein_unk_PA_PC_OL"
            water_ion_alch_ion_index = "Water_and_ions"

            for phase in ['coul', 'rest', 'vdw']:
                create_directories_with_MDP_pertuball_smaller_dumps(
                    phase, protein_membrane_ligand_index, water_ion_alch_ion_index
                )

            create_simulations_list()
            simulation_list_string = ' '.join([
                f"{phase}.{i:02d}/$STEP/" for phase in ['rest', 'coul', 'vdw']
                for i in range(21 if phase == 'vdw' else 12 if phase == 'rest' else 11)
            ])

            gen_hrex_submission_script(
                template_path=files("abfe.utils.cluster.job_script_template.abfe") / "job_complex_archer_hrex_template.sh",
                job_name=ligname,
                archer_nodes=self.archer_nodes,
                simulation_list_string=simulation_list_string
            )

            gen_hrex_submission_script(   
                template_path=files("abfe.utils.cluster.job_script_template.abfe") / "job_complex_archer_hrex_cont_template.sh",
                new_script_name='job_complex_archer_contd.sh',
                job_name=ligname,
                archer_nodes=self.archer_nodes,
                simulation_list_string=simulation_list_string
            )

            (complex_path / 'plumed.dat').touch()


            logging.info(f"ABFE folder setup completed for {ligname} replicate {suffix}")

        except Exception as e:
            logging.error(f"Error setting up ABFE folder for {ligand_dir.name} replicate {suffix}: {e}")
            raise


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Setup ABFE folders for multiple ligands and replicates.")
    parser.add_argument('--base_path', required=True, help='Base path containing ligand folders.')
    parser.add_argument('--ligands', nargs='+', required=True, help='Ligand folder names to process.')
    parser.add_argument('--num_replicates', type=int, default=3, help='Number of replicates to setup per ligand.')
    parser.add_argument('--archer_nodes', type=int, default=22, help='Number of nodes to use in scripts.')

    args = parser.parse_args()

    setup = ABFESetup(
        base_path=args.base_path,
        ligands=args.ligands,
        num_replicates=args.num_replicates,
        archer_nodes=args.archer_nodes
    )
    setup.setup()


if __name__ == "__main__":
    main()
