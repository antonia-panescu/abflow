import MDAnalysis as mda
from MDAnalysis.coordinates.GRO import GROWriter

def merge_gro_files(file1, file2, output_suffix='_crystal_water.gro'):
    u1 = mda.Universe(file1)
    u2 = mda.Universe(file2)
    combined = mda.Merge(u1.atoms, u2.atoms)
    combined.dimensions = u1.dimensions
    output_file = file1.replace('.gro', output_suffix)
    with GROWriter(output_file, n_atoms=combined.atoms.n_atoms) as W:
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

def addtoppar2top(filename="topol.top"):
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines = [l.replace('#include "posre.itp"', '#include "toppar/posre.itp"') for l in lines]
    with open(filename, 'w') as f:
        f.writelines(lines)

def copy_jobscripts_vanilla(template_dir, name, slots):
    for file_name in os.listdir(template_dir):
        path = os.path.join(template_dir, file_name)
        with open(path) as f:
            content = f.read()
        content = content.replace('$NAME$', name).replace('$NO_OF_NODES$', str(slots))
        with open(file_name, 'w') as f:
            f.write(content)
