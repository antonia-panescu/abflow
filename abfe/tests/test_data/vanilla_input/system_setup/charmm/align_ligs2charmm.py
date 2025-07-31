import os

def align_ligand_to_charmm(template_pml_path, ligand_sdf_path):
    """
    Replaces placeholders in a template PML file with the ligand SDF path and name,
    saves the tailored PML script, and executes it in the shell using PyMOL.
    
    Args:
        template_pml_path (str): Path to the template PML file.
        ligand_sdf_path (str): Path to the ligand SDF file.
    """
    # Extract the ligand SDF filename and base name
    ligand_file_name = os.path.basename(ligand_sdf_path)
    ligand_base_name = os.path.splitext(ligand_file_name)[0]
    
    # Read the template PML file
    with open(template_pml_path, "r") as template_file:
        pml_content = template_file.read()
    
    # Replace placeholders with the actual SDF path and name
    tailored_pml_content = pml_content.replace("$LIGAND_FILE_PATH", ligand_sdf_path)
    tailored_pml_content = tailored_pml_content.replace("$LIGAND_FILE_NAME", ligand_base_name+'.sdf')
    
    # Define the output tailored PML file path
    tailored_pml_path = "tailored_script.pml"
    
    # Write the tailored PML file
    with open(tailored_pml_path, "w") as tailored_file:
        tailored_file.write(tailored_pml_content)
    
    # Execute the tailored PML script using PyMOL in the shell
    command = f"~/pymol/installation/bin/pymol -cq {tailored_pml_path}"
    os.system(command)
    print(f"Executed tailored PML script: {tailored_pml_path}")

# Example usage
if __name__ == "__main__":
    # Replace with actual paths
    template_pml_path = "/biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/charmm/align_lig2charmm_template.pml"
    ligand_sdf_location = "/biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/Guo_et_al_SI/input_files/A2A_2020"
    
    # Run the function for each ligand SDF file in ligand sdf location
    for ligand_sdf in os.listdir(ligand_sdf_location):
        if not ligand_sdf.endswith(".sdf"):
            continue

        ligand_sdf_path = os.path.join(ligand_sdf_location, ligand_sdf)
        print(f"Aligning ligand SDF file: {ligand_sdf_path}")
        align_ligand_to_charmm(template_pml_path, ligand_sdf_path)
