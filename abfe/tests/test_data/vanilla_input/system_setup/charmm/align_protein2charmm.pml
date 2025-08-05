# Load the first PDB file
load /biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/protein_prep/Guo_A2A_2020_protein_w_waters.pdb, pdb1_object

# Load the second PDB file
load /biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/charmm/charmm_protein_only_step6.pdb, pdb2_object

# Align the merged object to the second PDB
align pdb1_object, pdb2_object

# Save the protein from the aligned merged object
select protein_only, pdb1_object and polymer.protein
save /biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/protein_prep/protein_aligned2charmm.pdb, protein_only

# Save the protein from the aligned merged object
select waters, pdb1_object and resn HOH
save /biggin/b211/reub0138/Projects/orexin/deflorian_A2A_v1/system_setup/protein_prep/crystal_waters_aligned2charmm.pdb, waters