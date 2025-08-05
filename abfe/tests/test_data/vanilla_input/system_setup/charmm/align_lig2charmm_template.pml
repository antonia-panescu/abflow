# Load the first PDB file
load /biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/Guo_et_al_SI/A2A_2020_protein_from_npt.pdb, pdb1_object

# Load the SDF file
load $LIGAND_FILE_PATH, sdf_object

# Merge the PDB1 and SDF into one object
create merged_object, pdb1_object or sdf_object

# Load the second PDB file
load /biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/charmm/charmm_protein_only_step6.pdb, pdb2_object

# Align the merged object to the second PDB
align merged_object, pdb2_object

# Save the 'resn UNK' from the aligned merged object
select resn_unk, merged_object and resn UNK
save ligands_aligned2charmm/$LIGAND_FILE_NAME, resn_unk
