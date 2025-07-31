import MDAnalysis as mda
from importlib.abc import Traversable
from MDAnalysis.coordinates.GRO import GROWriter
from pathlib import Path
from collections import Counter

def merge_gro_files(file1, file2, output_suffix='_crystal_water.gro'):
    u1 = mda.Universe(file1)
    u2 = mda.Universe(file2)
    combined = mda.Merge(u1.atoms, u2.atoms)
    combined.dimensions = u1.dimensions
    output_file = Path(file1).with_name(Path(file1).stem + output_suffix)
    with GROWriter(str(output_file), n_atoms=combined.atoms.n_atoms) as W:
        W.write(combined.atoms)
    print(f"Combined GRO file saved as {output_file}")

def add_water2topology(n=0):
    with open("topol.top", "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('SOL'):
            current_water_count = int(line.split()[1])
            print(f"Current number of water molecules: {current_water_count}")
            new_water_count = current_water_count + n
            print(f"New number of water molecules: {new_water_count}")
            lines[i] = f'SOL        {new_water_count}\n'
            break
    with open("topol.top", "w") as f:
        f.writelines(lines)



#def rebuild_molecule_counts_from_gro(gro_file: str) -> str:
    """
    Parses a GRO file using MDAnalysis and constructs a [ molecules ] block
    based on residue *segments* or custom logic.

    This version avoids listing individual amino acids (e.g., MET, GLY)
    and instead assumes the protein is treated as one moleculetype (e.g., 'Protein').

    Returns:
        str: Corrected [ molecules ] block.
    """
#   u = mda.Universe(gro_file)

#    mol_counts = Counter()

#    for res in u.residues:
#        resname = res.resname.strip()

#        if resname in {"SOL", "CL", "NA"}:
#            mol_counts[resname] += 1
#        elif res.segid.strip():  # Use segid if available
#            mol_counts[res.segid.strip()] += 1
#        elif res.resname == "unk":
#            mol_counts["UNK"] += 1
#        elif res.resid == 1:  # assume the first residue belongs to the protein
#            mol_counts["Protein"] += 1
#        else:
#            continue  # ignore individual residues like MET, GLY, etc.

#    print("🔬 Molecule counts (filtered):")
#    for k, v in mol_counts.items():
#        print(f"{k:<10}: {v}")

#    output = "[ molecules ]\n; Compound        #mols\n"
#    for mol, count in mol_counts.items():
#        output += f"{mol:<15}{count}\n"

#    return output




def addtoppar2top(filename="topol.top"):
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines = [l.replace('#include "posre.itp"', '#include "toppar/posre.itp"') for l in lines]
    with open(filename, 'w') as f:
        f.writelines(lines)


def copy_jobscripts_vanilla(template_dir: Traversable, name: str, slots: int):
    for template_file in template_dir.iterdir():
        if template_file.name.endswith('.sh'):
            with template_file.open('r') as f:
                content = f.read()

            content = content.replace('$NAME$', name).replace('$NO_OF_NODES$', str(slots))

            with open(template_file.name, 'w') as f_out:
                f_out.write(content)



#def update_topol_top_with_molecules_block(gro_file: str, top_file: str = "topol.top"):
#    mol_block = rebuild_molecule_counts_from_gro(gro_file)
#    with open(top_file, "r") as f:
#        lines = f.readlines()

#    with open(top_file, "w") as f:
#        inside_mol_block = False
#        for line in lines:
#            if line.strip().startswith("[ molecules ]"):
#                inside_mol_block = True
#                f.write(mol_block)
#            elif inside_mol_block:
#                if line.strip().startswith("[") and not line.strip().startswith("[ molecules ]"):
#                    inside_mol_block = False
#                    f.write(line)
#            elif not inside_mol_block:
#                f.write(line)
