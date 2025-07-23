import os
import shutil
import subprocess
import logging
from pathlib import Path
from rdkit import Chem


def run_gmx_command(command):      ## REFACTORED 
    """Run a GROMACS (gmx) command."""
    try:
        logging.info(f"Executing: {command}")
        subprocess.run(command, shell=True, check=True)
        logging.info("Success.")
    except subprocess.CalledProcessError as e:
        logging.error(f"GROMACS command failed: {e}")
        raise RuntimeError(f"GROMACS command failed: {e}")

def check_boresch(vanilla_folder: Path):        # REFACTORED
    """Check that vanilla folder and required files exist."""
    tpr_file = vanilla_folder / 'prod.tpr'
    xtc_file = vanilla_folder / 'prod.xtc'

    if not vanilla_folder.exists():
        raise FileNotFoundError(f"Missing folder: {vanilla_folder}")
    if not tpr_file.exists() or not xtc_file.exists():
        raise FileNotFoundError(f"Missing prod.tpr or prod.xtc in {vanilla_folder}")

def generate_boresch_restraints(vanilla_folder_name='vanilla', resname='unk', input_gro = 'npt.gro'):     
    """Creates Boresch restraints using the Boresch Restraint Generator"""
    boresch_folder = 'boresch'
    if not os.path.exists(boresch_folder):
        os.mkdir(boresch_folder)
        logging.info('Boresch folder created!')
    os.chdir(boresch_folder)
    
    check_boresch(vanilla_folder_name=vanilla_folder_name)
    
    logging.info('Creating a GROMACS tpr file...')
    #run_gmx_command(f'gmx grompp -f ../../{vanilla_folder_name}/prod.mdp -c ../../{vanilla_folder_name}/{input_gro} -r ../../{vanilla_folder_name}/{input_gro} -p ../../{vanilla_folder_name}/topol.top -o prod.tpr -n ../../{vanilla_folder_name}/index.ndx')
    run_gmx_command(f'gmx grompp -f ../../{vanilla_folder_name}/prod.mdp -c ../../{vanilla_folder_name}/{input_gro} -r ../../{vanilla_folder_name}/{input_gro} -p ../../{vanilla_folder_name}/topol_water_removed.top -o prod.tpr -n ../../{vanilla_folder_name}/index.ndx')


    logging.info('Running Boresch Restraint Generator with fit_fc...')
    script_path = os.path.abspath(os.path.expanduser('/biggin/b230/magd5710/Documents/Nithish_FEPA/Util/MDresgen/MDRestraintsGenerator/scripts/BoreschRestraintGMX.py'))
    command = [
        'python',
        script_path,
        '--top', 'prod.tpr',
        '--traj', f'../../{vanilla_folder_name}/prod.xtc',
        '--host_selection', 'protein and name CA',
        '--ligand_selection', f'resname {resname}'
    ]

    # Specify the output files
    stdout_file = "resgen_output.txt"
    stderr_file = "resgen_error.txt"

    # Run the command and save output and error to files
    try:
        with open(stdout_file, "w") as out_f, open(stderr_file, "w") as err_f:
            subprocess.run(command, stdout=out_f, stderr=err_f, check=True)
        logging.info(f"Output saved to: {stdout_file}")
        logging.info(f"Error saved to: {stderr_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running Boresch Restraint Generator: {e}")
        sys.exit(1)
    
    os.chdir('..')

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

def copy_complex(vanilla_folder_name='vanilla', copy_ff = True):
    """Copy necessary files to the complex folder and update topology."""
    # Create a folder called complex
    os.mkdir('complex')
    os.chdir('complex')
    # Copy the file 'ClosestRestraintFrame.gro','BoreschRestraint.top' from the boresch folder to the abfe/complex folder into the complex folder
    shutil.copy('../boresch/ClosestRestraintFrame.gro', 'complex.gro')
    shutil.copy('../boresch/BoreschRestraint.top', 'BoreschRestraint.top')
    #shutil.copy(f'../../{vanilla_folder_name}/topol.top', 'complex.top')
    shutil.copy(f'../../{vanilla_folder_name}/topol_water_removed.top', 'complex.top')
    shutil.copytree(f'../../{vanilla_folder_name}/toppar', 'toppar')
    if copy_ff:
        shutil.copytree('/biggin/b211/reub0138/Util/forcefields/amber99sb-star-ildn-mut.ff', 'amber99sb-star-ildn-mut.ff')
    add_boresch_to_topol('BoreschRestraint.top','complex.top')
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

def create_directories_with_MDP_pertuball(leg, pml_index, wiai_index):
    """Create directories and MDP files for different stages of simulation."""
    leg_windows = {
        "coul": (11, 12),
        "rest": (12, 0),
        "vdw": (21, 23)
    }

    windows, state_id_start = leg_windows.get(leg, (0, 0))
    state_id = state_id_start

    for i in range(windows):
        i_str = f"{i:02d}"
        directory_path = f"{leg}.{i_str}"
        os.makedirs(directory_path, exist_ok=True)
        
        for stage in ["enmin", "npt_b", "npt_pr", "nvt", "prod"]:
            stage_directory = os.path.join(directory_path, stage)
            os.makedirs(stage_directory, exist_ok=True)
            
            template_file = f"/biggin/b211/reub0138/Util/super_mdp_files_pertuball/{leg}.{stage}.super.mdp"
            output_file = f"{directory_path}/{stage}/{leg}.{stage}.{i_str}.mdp"

            with open(template_file, 'r') as f:
                template_content = f.read()
            
            modified_content = template_content.replace("Protein_LIG_POPC", pml_index)
            modified_content = modified_content.replace("Water_and_ions_B_CL", wiai_index)
            modified_content = modified_content.replace("<state>", str(state_id))

            with open(output_file, 'w') as f:
                f.write(modified_content)

        state_id += 1

def create_directories_with_MDP_pertuball_smaller_dumps(leg, pml_index, wiai_index):
    """Create directories and MDP files for different stages of simulation."""
    leg_windows = {
        "coul": (11, 12),
        "rest": (12, 0),
        "vdw": (21, 23)
    }

    windows, state_id_start = leg_windows.get(leg, (0, 0))
    state_id = state_id_start

    for i in range(windows):
        i_str = f"{i:02d}"
        directory_path = f"{leg}.{i_str}"
        os.makedirs(directory_path, exist_ok=True)
        
        for stage in ["enmin", "npt_b", "npt_pr", "nvt", "prod"]:
            stage_directory = os.path.join(directory_path, stage)
            os.makedirs(stage_directory, exist_ok=True)
            
            template_file = f"/biggin/b211/reub0138/Util/super_mdp_files_pertuball_smaller_dumps/{leg}.{stage}.super.mdp"
            output_file = f"{directory_path}/{stage}/{leg}.{stage}.{i_str}.mdp"

            with open(template_file, 'r') as f:
                template_content = f.read()
            
            modified_content = template_content.replace("Protein_LIG_POPC", pml_index)
            modified_content = modified_content.replace("Water_and_ions_B_CL", wiai_index)
            modified_content = modified_content.replace("<state>", str(state_id))

            with open(output_file, 'w') as f:
                f.write(modified_content)

        state_id += 1

def create_directories_with_soluble_MDP(leg, pml_index, wiai_index):
    """Create directories and MDP files for soluble simulations."""
    leg_windows = {
        "coul": 11,
        "rest": 12,
        "vdw": 21
    }

    windows = leg_windows.get(leg, 0)

    for i in range(windows):
        i_str = f"{i:02d}"
        directory_path = f"{leg}.{i_str}"
        os.makedirs(directory_path, exist_ok=True)
        
        for stage in ["enmin", "npt_b", "npt_pr", "nvt", "prod"]:
            stage_directory = os.path.join(directory_path, stage)
            os.makedirs(stage_directory, exist_ok=True)
            
            template_file = f"/biggin/b211/reub0138/Util/super_mdp_files_soluble/{leg}.{stage}.super.mdp"
            output_file = f"{directory_path}/{stage}/{leg}.{stage}.{i_str}.mdp"

            with open(template_file, 'r') as f:
                template_content = f.read()
            
            modified_content = template_content.replace("Protein_LIG_POPC", pml_index)
            modified_content = modified_content.replace("Water_and_ions_B_CL", wiai_index)
            modified_content = modified_content.replace("<state>", i_str)

            with open(output_file, 'w') as f:
                f.write(modified_content)


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

def gen_hrex_submission_script(template_path, job_name, archer_nodes, simulation_list_string, new_script_name = 'job_complex_archer.sh', qos='standard'):
    """Generate an Archer submission script from a template."""
    try:
        # Copy the template script to the current directory
        shutil.copy(template_path, new_script_name)

        # Calculate the Archer nodes setup
        tasks_per_node = 44 // archer_nodes
        cpus_per_task = 128 // tasks_per_node

        # Check for uneven division
        if 44 % archer_nodes != 0:
            logging.warning(f"Number of nodes ({archer_nodes}) does not divide 44 CPUs evenly.")
        if 128 % tasks_per_node != 0:
            logging.warning(f"Tasks per node ({tasks_per_node}) does not divide 128 CPUs per node evenly.")

        # Read the template script
        with open('job_complex_archer.sh', 'r') as file:
            filedata = file.read()

        # Replace placeholders with actual values
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

        # Write the modified script back to the file
        with open(new_script_name, 'w') as file:
            file.write(filedata)

        logging.info("Archer submission script generated successfully.")
    except Exception as e:
        logging.error(f"Error generating Archer submission script: {e}")
        sys.exit(1)