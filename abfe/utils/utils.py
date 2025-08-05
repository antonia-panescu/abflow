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

#def add_water2topology(n=0):
#    with open("topol.top", "r") as f:
#        lines = f.readlines()
#    for i, line in enumerate(lines):
#        if line.startswith('SOL'):
#            current_water_count = int(line.split()[1])
#            print(f"Current number of water molecules: {current_water_count}")
#            new_water_count = current_water_count + n
#            print(f"New number of water molecules: {new_water_count}")
#            lines[i] = f'SOL        {new_water_count}\n'
#            break
#    with open("topol.top", "w") as f:
#        f.writelines(lines)



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




def add_water2topology(n: int = 0):
    """
    Increment the water molecule count in topol.top by n.
    Recognises common solvent residue names (SOL, TP3, WAT, HOH).
    """
    with open("topol.top", "r") as f:
        lines = f.readlines()

    solvent_keys = ("SOL", "TP3", "WAT", "HOH")
    for i, line in enumerate(lines):
        tokens = line.split()
        if tokens and tokens[0] in solvent_keys:
            current_water_count = int(tokens[1])
            new_count = current_water_count + n
            # Preserve the solvent name and column width
            lines[i] = f"{tokens[0]:<10}{new_count}\n"
            break

    with open("topol.top", "w") as f:
        f.writelines(lines)




def rebuild_molecule_counts_from_gro(gro_file: str) -> str:
    """
    Parse a .gro file and build a new [ molecules ] block by counting residues.
    Water names (SOL, TP3, WAT, HOH) are all treated as solvent.
    Residue segids are used to group membrane lipids.
    """
    u = mda.Universe(gro_file)
    mol_counts = Counter()

    for res in u.residues:
        name = res.resname.strip()
        if name in {"SOL", "TP3", "WAT", "HOH"}:
            mol_counts[name] += 1
        elif res.segid.strip():
            mol_counts[res.segid.strip()] += 1
        elif res.resname == "unk":
            mol_counts["UNK"] += 1
        elif res.resid == 1:
            mol_counts["Protein"] += 1

    output = "[ molecules ]\n; Compound        #mols\n"
    for mol, count in mol_counts.items():
        output += f"{mol:<15}{count}\n"
    return output

def update_topol_top_with_molecules_block(gro_file: str, top_file: str = "topol.top"):
    """
    Replace the [ molecules ] section of top_file with counts derived from gro_file.
    """
    mol_block = rebuild_molecule_counts_from_gro(gro_file)
    with open(top_file, "r") as f:
        lines = f.readlines()

    with open(top_file, "w") as f:
        inside_mol_block = False
        for line in lines:
            stripped = line.strip()
            # When we encounter the start of the old [ molecules ] block,
            # write the rebuilt block instead.
            if stripped.startswith("[ molecules ]"):
                inside_mol_block = True
                f.write(mol_block)
                continue
            # Skip lines until we leave the old [ molecules ] block.
            if inside_mol_block:
                if stripped.startswith("[") and not stripped.startswith("[ molecules ]"):
                    inside_mol_block = False
                    f.write(line)
                # Otherwise continue skipping.
                continue
            f.write(line)

import os
from contextlib import contextmanager
from pathlib import Path

@contextmanager
def working_directory(path: Path):
    """
    Temporarily change the working directory inside a `with` block.
    Restores the original working directory upon exit.
    """
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
