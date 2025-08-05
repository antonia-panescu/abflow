import os
import shutil
import subprocess
import sys
import logging
from pathlib import Path
from rdkit import Chem
from importlib.abc import Traversable
import importlib.resources as pkg_resources
import numpy as _np
from itertools import chain


import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from MDRestraintsGenerator.search import find_ligand_atoms, FindHostAtoms
from MDRestraintsGenerator.restraints import FindBoreschRestraint


#sys.path.insert(0, '/biggin/b230/magd5710/Documents/Nithish_FEPA/Util/MDresgen/MDRestraintsGenerator')


from abfe.utils import data 
from abfe.utils.data import forcefield_dir 
from abfe.utils.data.mdps import super_mdp_files_pertuball, super_mdp_files_pertuball_smaller_dumps, super_mdp_files_soluble 



def run_gmx_command(command):      ## REFACTORED 
    """Run a GROMACS (gmx) command."""
    try:
        logging.info(f"Executing: {command}")
        subprocess.run(command, shell=True, check=True)
        logging.info("Success.")
    except subprocess.CalledProcessError as e:
        logging.error(f"GROMACS command failed: {e}")
        raise RuntimeError(f"GROMACS command failed: {e}")


def check_boresch(vanilla_folder):  # Accepts str or Path
    """Check that vanilla folder and required files exist."""
    vanilla_folder = Path(vanilla_folder)  # Ensure it's a Path object

    tpr_file = vanilla_folder / 'prod.tpr'
    xtc_file = vanilla_folder / 'prod.xtc'

    if not vanilla_folder.exists():
        raise FileNotFoundError(f"Missing folder: {vanilla_folder}")
    if not tpr_file.exists() or not xtc_file.exists():
        raise FileNotFoundError(f"Missing prod.tpr or prod.xtc in {vanilla_folder}")

def _flatten_indices(nested):
    """Recursively flatten lists/tuples/ndarrays into a flat list of ints."""
    flat = []
    if isinstance(nested, (list, tuple, _np.ndarray)):
        for x in nested:
            flat.extend(_flatten_indices(x))
    else:
        flat.append(int(nested))
    return flat


def generate_boresch_restraints(
    vanilla_folder_name: str = "vanilla",
    resname:           str = "unk",
    input_gro:         str = "npt.gro",
    gmx_exec:         str = "gmx"
):
    """
    Create Boresch restraints via the MDRestraintsGenerator API.
    Results go into <vanilla_folder>/boresch/.
    """
    vanilla = Path(vanilla_folder_name).resolve()
    boresch = vanilla / "boresch"
    boresch.mkdir(exist_ok=True)

    tpr = boresch / "prod.tpr"
    xtc = vanilla / "prod.xtc"

    # 1) build the tpr
    cmd = (
        f"{gmx_exec} grompp "
        f"-f {vanilla/'prod.mdp'} "
        f"-c {vanilla/input_gro} "
        f"-r {vanilla/input_gro} "
        f"-p {vanilla/'topol_water_removed.top'} "
        f"-n {vanilla/'index.ndx'} "
        f"-o {tpr}"
    )
    logging.info("Creating GROMACS tpr file…")
    subprocess.run(cmd, shell=True, check=True)

    # 2) call the API
    logging.info("Running Boresch restraint generation via API…")
    try:
        u = mda.Universe(str(tpr), str(xtc))

        # sanity‐check ligand residue name
        residues = {r.resname for r in u.residues}
        if resname not in residues:
            raise ValueError(f"resname '{resname}' not in universe; found {sorted(residues)}")

        # pick ligand atoms → raw indices or AtomGroup
        raw_lig = find_ligand_atoms(u, f"resname {resname}")
        if isinstance(raw_lig, AtomGroup):
            lig = raw_lig
            lig_idxs = lig.indices.tolist()
        else:
            lig_idxs = _flatten_indices(raw_lig)
            lig = u.atoms[lig_idxs]

        # pick host atoms by index list (must be ints, not AtomGroup)
        raw_host = FindHostAtoms(u, lig_idxs, "protein and name CA").run()
        if isinstance(raw_host, AtomGroup):
            host = raw_host
        else:
            host_idxs = _flatten_indices(raw_host)
            host = u.atoms[host_idxs]

        # build & write the restraint
        bo = FindBoreschRestraint(u, lig, host)
        bo.run()
        bo.write(prefix="BoreschRestraint", directory=str(boresch))

        logging.info("Boresch restraints successfully generated.")
    except Exception:
        logging.exception("Error generating Boresch restraints")
        sys.exit(1)


def check_line(filename, text):
    """Check if a specific text is present in a file."""
    if not os.path.isfile(filename):
        logging.error(f"File '{filename}' not found.")
        return False

    with open(filename, 'r') as file:
        for line in file:
            if text in line:
                return True
    return False
    
def add_boresch_to_topol(source_file, destination_file):
    """Add Boresch restraints from source_file to destination_file if not already present."""
    if check_line(destination_file, '; restraints'):
        logging.info("The 'topol.top' file already contains the line '; restraints'.")
    else:
        logging.info("The 'topol.top' file does not contain the line '; restraints'. Adding...")
        try:
            with open(source_file, 'r') as source:
                source_contents = source.read()
            with open(destination_file, 'a') as destination:
                destination.write('\n; restraints\n')
                destination.write(source_contents)
            logging.info("Boresch restraints added successfully.")
        except IOError as e:
            logging.error(f"Error reading or writing files: {e}")
            sys.exit(1)


def copy_complex(vanilla_folder_name='vanilla', copy_ff=True):
    """
    Copy necessary files to the 'complex' folder and update topology.
    
    Parameters:
        vanilla_folder_name (str or Path): Folder containing vanilla system.
        copy_ff (bool): Whether to copy the forcefield directory.
    """
    vanilla_folder = Path(vanilla_folder_name)
    complex_dir = Path('complex')
    complex_dir.mkdir(exist_ok=True)

    # Copy necessary files
    shutil.copy('boresch/ClosestRestraintFrame.gro', complex_dir / 'complex.gro')
    shutil.copy('boresch/BoreschRestraint.top', complex_dir / 'BoreschRestraint.top')
    shutil.copy(vanilla_folder / 'topol_water_removed.top', complex_dir / 'complex.top')
    shutil.copytree(vanilla_folder / 'toppar', complex_dir / 'toppar', dirs_exist_ok=True)

    # Copy forcefield from package data
    if copy_ff:
        ff_src = pkg_resources.files(forcefield_dir) / 'amber99sb-star-ildn-mut.ff'
        ff_dst = complex_dir / 'amber99sb-star-ildn-mut.ff'
        shutil.copytree(ff_src, ff_dst, dirs_exist_ok=True)

    # Update topology file
    add_boresch_to_topol(complex_dir / 'BoreschRestraint.top', complex_dir / 'complex.top')
    logging.info("Files copied and topology updated successfully.")


def calculate_total_charge(sdf_file):
    """Calculate the total formal charge of a molecule from an SDF file."""
    # Load the molecule from the SDF file
    supplier = Chem.SDMolSupplier(sdf_file)
    
    # Check if there is a molecule in the SDF file
    if not supplier or len(supplier) == 0:
        raise ValueError("No molecules found in the SDF file.")
    
    mol = supplier[0]  # Assuming there is only one molecule in the file
    
    # Calculate the total charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    logging.info(f'Ligand charge: {total_charge}')
    
    return total_charge

def create_fep_system(lig_name: str, vanilla_folder: Path, force=False):        ## REFACTORED
    """Create a Free Energy Perturbation (FEP) system for a given ligand."""
    try:
        # Calculate the ligand charge
        lig_sdf = vanilla_folder / 'smirnoff' / f"{lig_name}.sdf"
        lig_ch = calculate_total_charge(lig_sdf)
        if lig_ch != 0:
            logging.warning(f"Ligand {lig_name} has a non-zero charge of {lig_ch} e.")
            if not force:
                user_input = input(f"Ligand {lig_name} charge = {lig_ch} e. Proceed? (y/n): ")
                if user_input.lower() != 'y':
                    logging.info("Aborting FEP system creation by user choice.")
                    return

        # Define the command parameters
        alch_mthd = 'co-annihilation'
        lig_resname = "unk"
        lig_itp = f"toppar/{lig_name}.itp"
        lig_atomtype_itp = f"toppar/{lig_name}_atomtypes.itp"
        script_path = Path(__file__).parent / "add_alchemical_ion_v3_pertuball.py"
        command = [
            'python',
            script_path,
            '--lig_charge', lig_ch,
            '--alch_ion_method', alch_mthd,
            '--lig', lig_resname,
            '--lig_itp', lig_itp,
            '--atomtype_itp', lig_atomtype_itp
        ]

        # Run the command
        subprocess.run(command, check=True)
        print('Ligand itp:', lig_itp)
        print('Ligand atomtype itp:', lig_atomtype_itp)
    except subprocess.CalledProcessError as e:
        print(f"ERROR with Alchemical Ion v3: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)

def create_index(input_file, output_file):
    """
    Create an index file using gmx make_ndx command.
    
    Parameters:
        input_file (str): Input structure file (.gro, .pdb, etc.).
        output_file (str): Output index file.
    """
    try:
        # Run the gmx make_ndx command
        logging.info(f"Creating index file: create Protein_lig_mem and water_ion_alchion groups")
        command = f'echo "1|13|23|24|25\n q\n" | gmx make_ndx -f {input_file} -o {output_file}'
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Index file {output_file} created successfully.")
        
        # Copy index.ndx to index_coul.ndx, index_vdw.ndx, and index_rest.ndx
        for suffix in ['coul', 'vdw', 'rest']:
            shutil.copy(output_file, f'index_{suffix}.ndx')
        logging.info("Index files copied to index_coul.ndx, index_vdw.ndx, and index_rest.ndx successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while creating index file: {e}")
        sys.exit(1)
    except IOError as e:
        logging.error(f"Error copying index files: {e}")
        sys.exit(1)

def create_directories_with_MDP_pertuball(leg: str, pml_index: str, wiai_index: str):
    """
    Create FEP simulation directories and fill them with modified .mdp files
    from the packaged super_mdp_files_pertuball directory.

    Parameters:
        leg (str): Perturbation leg ('coul', 'rest', 'vdw').
        pml_index (str): Index group for protein/membrane/ligand.
        wiai_index (str): Index group for water/ions/alchemical ion.
    """
    leg_windows = {
        "coul": (11, 12),
        "rest": (12, 0),
        "vdw": (21, 23)
    }

    windows, state_id_start = leg_windows.get(leg, (0, 0))
    state_id = state_id_start
    stages = ["enmin", "npt_b", "npt_pr", "nvt", "prod"]

    for i in range(windows):
        i_str = f"{i:02d}"
        dir_base = Path(f"{leg}.{i_str}")
        dir_base.mkdir(exist_ok=True)

        for stage in stages:
            stage_dir = dir_base / stage
            stage_dir.mkdir(exist_ok=True)

            template_filename = f"{leg}.{stage}.super.mdp"
            output_filename = stage_dir / f"{leg}.{stage}.{i_str}.mdp"

            try:
                with pkg_resources.files(super_mdp_files_pertuball).joinpath(template_filename).open("r") as f:
                    template_content = f.read()

                modified_content = (
                    template_content
                    .replace("Protein_LIG_POPC", pml_index)
                    .replace("Water_and_ions_B_CL", wiai_index)
                    .replace("<state>", str(state_id))
                )

                with open(output_filename, "w") as f_out:
                    f_out.write(modified_content)

                logging.info(f"Written: {output_filename}")
            except FileNotFoundError:
                logging.error(f"Template file not found: {template_filename}")
                continue

        state_id += 1

def create_directories_with_MDP_pertuball_smaller_dumps(leg: str, pml_index: str, wiai_index: str):
    """
    Create FEP simulation directories using MDPs with smaller dump frequencies,
    sourced from packaged `super_mdp_files_pertuball_smaller_dumps`.

    Parameters:
        leg (str): 'coul', 'rest', or 'vdw'.
        pml_index (str): Index group name for protein/membrane/ligand.
        wiai_index (str): Index group name for water/ions/alchemical ion.
    """
    leg_windows = {
        "coul": (11, 12),
        "rest": (12, 0),
        "vdw": (21, 23)
    }

    windows, state_id_start = leg_windows.get(leg, (0, 0))
    state_id = state_id_start
    stages = ["enmin", "npt_b", "npt_pr", "nvt", "prod"]

    for i in range(windows):
        i_str = f"{i:02d}"
        dir_base = Path(f"{leg}.{i_str}")
        dir_base.mkdir(exist_ok=True)

        for stage in stages:
            stage_dir = dir_base / stage
            stage_dir.mkdir(exist_ok=True)

            template_filename = f"{leg}.{stage}.super.mdp"
            output_filename = stage_dir / f"{leg}.{stage}.{i_str}.mdp"

            try:
                with pkg_resources.files(super_mdp_files_pertuball_smaller_dumps).joinpath(template_filename).open("r") as f:
                    template_content = f.read()

                modified_content = (
                    template_content
                    .replace("Protein_LIG_POPC", pml_index)
                    .replace("Water_and_ions_B_CL", wiai_index)
                    .replace("<state>", str(state_id))
                )

                with open(output_filename, "w") as f_out:
                    f_out.write(modified_content)

                logging.info(f"Written: {output_filename}")
            except FileNotFoundError:
                logging.error(f"Missing MDP template: {template_filename}")
                continue

        state_id += 1

def create_directories_with_soluble_MDP(leg: str, pml_index: str, wiai_index: str):
    """
    Create FEP simulation directories for soluble systems using MDP templates
    from the `super_mdp_files_soluble` package data.

    Parameters:
        leg (str): Perturbation leg ('coul', 'rest', 'vdw').
        pml_index (str): Index group name for protein/ligand.
        wiai_index (str): Index group name for water/ions/alchemical ion.
    """
    leg_windows = {
        "coul": 11,
        "rest": 12,
        "vdw": 21
    }

    windows = leg_windows.get(leg, 0)
    stages = ["enmin", "npt_b", "npt_pr", "nvt", "prod"]

    for i in range(windows):
        i_str = f"{i:02d}"
        dir_base = Path(f"{leg}.{i_str}")
        dir_base.mkdir(exist_ok=True)

        for stage in stages:
            stage_dir = dir_base / stage
            stage_dir.mkdir(exist_ok=True)

            template_filename = f"{leg}.{stage}.super.mdp"
            output_filename = stage_dir / f"{leg}.{stage}.{i_str}.mdp"

            try:
                with pkg_resources.files(super_mdp_files_soluble).joinpath(template_filename).open("r") as f:
                    template_content = f.read()

                modified_content = (
                    template_content
                    .replace("Protein_LIG_POPC", pml_index)
                    .replace("Water_and_ions_B_CL", wiai_index)
                    .replace("<state>", i_str)
                )

                with open(output_filename, "w") as f_out:
                    f_out.write(modified_content)

                logging.info(f"Wrote: {output_filename}")
            except FileNotFoundError:
                logging.error(f"Missing soluble MDP template: {template_filename}")
                continue


def create_simulations_list():
    """Create a simulations_list.txt file with simulation directories."""
    leg_windows = {
        "coul": 11,
        "rest": 12,
        "vdw": 21
    }

    with open("simulations_list.txt", "w") as file:
        for leg, windows in leg_windows.items():
            for i in range(windows):
                file.write(f"{leg}.{i:02d}\n")

def gen_hrex_submission_script(template_path: Traversable, job_name, archer_nodes, simulation_list_string, new_script_name='job_complex_archer.sh', qos='standard'):
    """Generate an Archer submission script from a packaged template."""
    try:
        # Calculate Archer layout
        tasks_per_node = 44 // archer_nodes
        cpus_per_task = 128 // tasks_per_node

        if 44 % archer_nodes != 0:
            logging.warning(f"Number of nodes ({archer_nodes}) does not divide 44 CPUs evenly.")
        if 128 % tasks_per_node != 0:
            logging.warning(f"Tasks per node ({tasks_per_node}) does not divide 128 CPUs per node evenly.")

        # Read the template content
        with template_path.open('r') as f:
            filedata = f.read()

        # Replace placeholders
        replacements = {
            '$JOBNAME': job_name,
            '$NO_OF_NODES': str(archer_nodes),
            '$TASKS_PER_NODE': str(tasks_per_node),
            '$CPUS_PER_TASK': str(cpus_per_task),
            '$WINDOW_LIST_STR': simulation_list_string.replace('/$STEP/', ''),
            '$MULTIDIR_LIST_STR_ENMIN': simulation_list_string.replace('$STEP', 'enmin'),
            '$MULTIDIR_LIST_STR_NVT': simulation_list_string.replace('$STEP', 'nvt'),
            '$MULTIDIR_LIST_STR_NPT_B': simulation_list_string.replace('$STEP', 'npt_b'),
            '$MULTIDIR_LIST_STR_NPT_PR': simulation_list_string.replace('$STEP', 'npt_pr'),
            '$MULTIDIR_LIST_STR_PROD': simulation_list_string.replace('$STEP', 'prod'),
            '$QOS': qos
        }

        for key, value in replacements.items():
            filedata = filedata.replace(key, value)

        # Write to output script file
        with open(new_script_name, 'w') as f_out:
            f_out.write(filedata)

        logging.info(f"Submission script '{new_script_name}' generated successfully.")

    except Exception as e:
        logging.error(f"Error generating Archer submission script: {e}")
        sys.exit(1)