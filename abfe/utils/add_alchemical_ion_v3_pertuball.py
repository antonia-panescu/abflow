'''
This is a simple, self-contained script to generate certain input files for ABFEs.
Although it will work with neutral ligands, it is specifically designed to work with 
charged ligands, when an alchemical ion is to be used.

Example usage:
            python add_alchemical_ion.py --lig_charge -1 --alch_ion_method charge-transfer

You can specify:
    The name and charge of the ligand
    The alchemical ion method to be used
    The decoupling method to be used

Expected inputs: 
    gro file
    lig.itp file
    lig_atomtype.itp
    topol.top 
Expected outputs: 
    top files for rest, coul, and vdw
    gro files for rest, coul, and vdw
    lig itp files for rest, coul, and vdw
    atomtype itp files for rest, coul, and vdw

'''

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import argparse
import warnings
import os
import shutil
from pprint import pprint
import pandas as pd
import math


if __name__ == "__main__":
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--lig', default='LIG',
                            help=('Name of the ligand'))
        parser.add_argument('--lig_charge', default=0, type=int,
                            help=('Charge of the ligand'))
        parser.add_argument('--alch_ion_method', default=None,
                            help=('method for maintaining charge neutrality'
                                '\nAcceptable values are:{acceptable_alch_ion_methods}'))
        parser.add_argument('--top', default='complex.top',
                            help='the input topology')
        parser.add_argument('--gro', default='complex.gro',
                            help='the input gro file')
        parser.add_argument('--atomtype_itp', default='lig_atomtype.itp',
                            help='the input atomtypes itp\nThis must be referenced in top file')
        parser.add_argument('--lig_itp', default='lig.itp',
                            help='the input ligand itp\nThis must be referenced in top file')
        parser.add_argument('--decouple_method', default='dummy-pertuball',
                            help=(f'the type of decoupling to do\nThis will change what type of itps are generated.'
                                '\nAcceptable values are:{acceptable_decouple_args}'))
        parser.add_argument('--gen_itps', default=True,
                            help='Whether to generate itps for ABFEs\nAcceptable options are: "True", "False"')
        args = parser.parse_args()
        return args

    def check_args(args):
        if args.lig_charge != 0:
            assert args.alch_ion_method in acceptable_alch_ion_methods, f'You are providing a charged ligand.\n"--alch_ion_method" must be one of {acceptable_alch_ion_methods}'
        elif args.lig_charge == 0 and args.alch_ion_method != None:
            args.alch_ion_method = None
            print('WARNING:\nThe ligand is neutral but you have requested a counterion.\nThis does not make sense, so we have set "alch_ion_method=None"')
        assert args.decouple_method in acceptable_decouple_args, f'Print "--decouple_method must be one of {acceptable_decouple_args}'
        if args.decouple_method == 'ci_no':
            print('WARNING:\nAre you sure you want to set couple-intramol=no?\nThis is known to cause issues in gromacs (as of gromacs 2023 beta)\n')
        assert os.path.isfile(args.top), f'No file with name {args.top}'
        assert os.path.isfile(args.gro), f'No file with name {args.gro}'

    def get_system_charge(gro, alch_ion, ion, ion_mass):
        u = mda.Universe(gro)
        charge_res = {'LYS': 1, 'ARG': 1, 'GLU': -1,
                      'NA': 1, 'CL': -1, 'K': 1,'HIP':1,'CLA':-1,'SOD':1,
                       args.lig: int(args.lig_charge)}
        
        if args.lig_charge < 0:
            charge_res[alch_ion] = charges['negative'][args.alch_ion_method]['rest']['a']
        if args.lig_charge > 0:
            charge_res[alch_ion] = int(
                charges['positive'][args.alch_ion_method]['rest']['a'])

        total_charge = 0
        for residue in u.residues:
            if residue.resname in charge_res.keys():
                print(f'Adding charge for residue {residue.resname} with charge {charge_res[residue.resname]}')
                total_charge += charge_res[residue.resname]
            if residue.resname == 'ASP' and len(residue.atoms.names) > 12:
                total_charge -= 0
                print(f'Adding charge for residue {residue.resname} with charge 0 {len(residue.atoms.names)}')
            elif residue.resname == 'ASP' and len(residue.atoms.names) <= 12:
                total_charge -= 1
                print(f'Adding charge for residue {residue.resname} with charge -1 {len(residue.atoms.names)}')
            if residue.resname == 'HIS' and len(residue.atoms.names) > 17:
                total_charge += 1
                print(f'Adding charge for residue {residue.resname} with charge 1 (HIS with 18 atoms)')
        total_charge = total_charge 
        print('Total charge of the system from complex.gro:', total_charge)
        return total_charge

    def create_charges():
        charges = {}
        for charge in ['negative', 'positive']:
            charges[charge] = {}
            for method in ['co-annihilation', 'charge-transfer']:
                charges[charge][method] = {}
                for stage in ['coul', 'rest', 'vdw']:
                    charges[charge][method][stage] = {}
                if charge == 'negative' and method == 'co-annihilation':
                    charges[charge][method]['rest']['a'] = 1
                    charges[charge][method]['rest']['b'] = 1
                    charges[charge][method]['coul']['a'] = 1
                    charges[charge][method]['coul']['b'] = 0
                    charges[charge][method]['vdw']['a'] = 0
                    charges[charge][method]['vdw']['b'] = 0
                if charge == 'positive' and method == 'co-annihilation':
                    charges[charge][method]['rest']['a'] = -1
                    charges[charge][method]['rest']['b'] = -1
                    charges[charge][method]['coul']['a'] = -1
                    charges[charge][method]['coul']['b'] = 0
                    charges[charge][method]['vdw']['a'] = 0
                    charges[charge][method]['vdw']['b'] = 0
                if charge == 'negative' and method == 'charge-transfer':
                    charges[charge][method]['rest']['a'] = 0
                    charges[charge][method]['rest']['b'] = 0
                    charges[charge][method]['coul']['a'] = 0
                    charges[charge][method]['coul']['b'] = -1
                    charges[charge][method]['vdw']['a'] = -1
                    charges[charge][method]['vdw']['b'] = -1
                if charge == 'positive' and method == 'charge-transfer':
                    charges[charge][method]['rest']['a'] = 0
                    charges[charge][method]['rest']['b'] = 0
                    charges[charge][method]['coul']['a'] = 0
                    charges[charge][method]['coul']['b'] = 1
                    charges[charge][method]['vdw']['a'] = 1
                    charges[charge][method]['vdw']['b'] = 1

        return charges

    def get_counter_ion_info(lig_charge, alch_ion_method):
        if lig_charge > 0 and alch_ion_method == 'co-annihilation':
            alch_ion = 'B_CL'
            ion = 'Cl'
            ion_mass = 35.45
        elif lig_charge < 0 and alch_ion_method == 'co-annihilation':
            alch_ion = 'B_NA'
            ion = 'Na'
            ion_mass = 22.99
        elif lig_charge < 0 and alch_ion_method == 'charge-transfer':
            alch_ion = 'B_CL'
            ion = 'Cl'
            ion_mass = 35.45
        elif lig_charge > 0 and alch_ion_method == 'charge-transfer':
            alch_ion = 'B_NA'
            ion = 'Na'
            ion_mass = 22.99
        else:
            alch_ion = None
            ion = None
            ion_mass = None

        return alch_ion, ion, ion_mass

    def write_file(file,list_of_lines):
        outfile = open(file,'w')
        for line in list_of_lines:
            if '\n' not in line:
                outfile.write(f'{line}\n')
            else:
                outfile.write(f'{line}')
        outfile.close()

    def add_counter_ion_to_gro(gro, stage, lig_charge, alch_ion_method, alch_ion, ion, ion_mass):
        if lig_charge == 0:
            shutil.copyfile(gro,f'{gro[:-4]}_{stage}.gro')
            return 0,0

        elif lig_charge > 0:
            charge_a = charges['positive'][alch_ion_method][stage]['a']
            charge_b = charges['positive'][alch_ion_method][stage]['b']
        elif lig_charge < 0:
            charge_a = charges['negative'][alch_ion_method][stage]['a']
            charge_b = charges['negative'][alch_ion_method][stage]['b']

        number_of_ions = abs(lig_charge)
        u = mda.Universe(gro)
        ions = u.select_atoms(f'resname {ion.upper()}')
        lig = u.select_atoms(f'resname {args.lig}')
        distances = list(distance_array(
            lig.center_of_geometry(), ions.positions)[0])

        for index, i in enumerate(range(number_of_ions)):
            positions = [[0, 0, 0], [u.dimensions[0]/2, 0, 0],
                         [0, u.dimensions[1]/2, 0]]
            ions[distances.index(max(distances))].residue.resname = alch_ion
            ions[distances.index(max(distances))].position = positions[index]
            distances[distances.index(max(distances))] = 0

        # we must remove an oppositely charged ion to remain neutral
        if alch_ion_method == 'charge-transfer' and lig_charge > 0:
            sacrificial = 'CL'
        elif alch_ion_method == 'charge-transfer' and lig_charge < 0:
            sacrificial = 'NA'
        if alch_ion_method == 'charge-transfer' and lig_charge != 0:
            ions = u.select_atoms(f'resname {sacrificial}')
            lig = u.select_atoms(f'resname {args.lig}')
            distances = list(distance_array(
                lig.center_of_geometry(), ions.positions)[0])
            for index, i in enumerate(range(number_of_ions)):
                ions[distances.index(max(distances))].residue.resname = 'REM'
                distances[distances.index(max(distances))] = 0

        new_u = mda.Merge(u.select_atoms(
            f'not resname {alch_ion} REM'), u.select_atoms(f'resname {alch_ion}'))
        new_u.dimensions = u.dimensions
        new_u.atoms.write(f'{gro[:-4]}_{stage}.gro')
    
        assert get_system_charge(f'{gro[:-4]}_{stage}.gro', alch_ion, ion, ion_mass) == 0, f'{gro[:-4]}_{stage}.gro is not neutral. This should not happen - my bad.'
        return charge_a, charge_b

    def add_counter_ion_to_top(alch_ion, ion, ion_mass, charge_a, charge_b, stage, gro, top):
        if args.lig_charge == 0:
            #print('Adding topology for neutral ligand:')
            out_lines = []
            for line in open(top).readlines():
                if args.atomtype_itp in line:
                    line = line.replace(
                        args.atomtype_itp, f'{args.atomtype_itp[:-4]}_{stage}.itp')
                if args.lig_itp in line:
                    out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+stage+'.itp'
                    line = line.replace(args.lig_itp,out_name)
                out_lines.append(line)
            #sys.exit()
            print(f'Writing file: {top[:-4]}_{stage}.top')
            write_file(f'{top[:-4]}_{stage}.top', out_lines)
        else:
            alch_ion_count = 0 # this is to check that the alchemical ion has been added
            alch_ion_top_string = (f'\n[ moleculetype ]\n{alch_ion} 1\n[ atoms ]\n1 {ion} 1 {ion.upper()} {ion.upper()}'
                f' 1 {charge_a} {ion_mass} {ion} {charge_b} {ion_mass}\n[ position_restraints ]\n 1 1 1000 1000 1000\n\n')
            u = mda.Universe(gro)
            ref_residues = ['SOL', 'NA', 'CL', alch_ion]
            top_pos = ''
            out_lines = []
            for index, line in enumerate(open(args.top).readlines()):
                if line[0] == '[':
                    top_pos = line[0:11]
                if args.atomtype_itp in line:
                    atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+stage+'.itp'
                    line = line.replace(
                        args.atomtype_itp, atom_out_name)
                #print(args.lig_itp)
                if args.lig_itp in line:
                    out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+stage+'.itp'
                    line = line.replace(args.lig_itp,out_name)
                if '[ system ]' in line:
                    out_lines.append(alch_ion_top_string)
                if top_pos == '[ molecules':
                    if len(line.split()) > 0 and line[0] != ';':
                        atomtype = line.split()[0]
                        if atomtype in ref_residues:
                            line = f'{atomtype} {len(u.select_atoms(f"resname {atomtype}").residues)}\n'
                    if line == '\n':
                        line = ''
                    if '; restraints' in line:
                        line = ''
                if '[ intermolecular_interactions ]' in line:
                    out_lines.append(f'{alch_ion} {len(u.select_atoms(f"resname {alch_ion}").residues)}\n; restraints\n')
                    alch_ion_count=1
                out_lines.append(line)
            if alch_ion_count == 0:
                print(f'WARNING: You seem to be missing the "[ intermolecular_interactions ]" section from your {args.top} file')
                out_lines.append(f'{alch_ion} {len(u.select_atoms(f"resname {alch_ion}").residues)}\n')
            write_file(f'{top[:-4]}_{stage}.top',out_lines)

    def gen_dummy_line(line, stage):
        line = line.split()
        line_fmt = ("%5d %10s %6d %6s %6s %6d %10.8f %10.6f %6s %10.8f %10.6f")
        if stage == 'vdw':
            atom_charge = 0.0
        elif stage == 'coul':
            #print(line)
            atom_charge = float(line[6])
        #print(stage)
        new_line = (line_fmt % (int(line[0]), line[1], int(line[2]),
                                line[3], line[4], int(line[5]),
                                atom_charge, float(line[7]), f"d{line[1]}",
                                0.0, float(line[7])))
        return new_line

    def gen_dummy_itps(lig_itp, stage, outfile):
        out_lines=[]
        section = ''
        for line in open(lig_itp).readlines():
            if line[0] == '[':
                section = line[:-1]
            if '[ atoms ]' in section:
                if line[0] != ';' and len(line) > 5 and '[ atoms ]' not in line:
                    line = gen_dummy_line(line, stage)
            out_lines.append(line)
        write_file(outfile,out_lines)
        return

    def gen_dummy_atomtypes(atomtype_itp, stage, outfile):
        colnames = ['names', 'atnum', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        df = pd.read_csv(atomtype_itp, sep="\s+", names=colnames,
                         comment=';', on_bad_lines='warn', skiprows=1)
        with open(outfile, 'w') as writer:
            writer.write("[ atomtypes ]\n")
            writer.write(
                "; name    at.num    mass      charge   ptype     sigma         epsilon\n")
            row_fmt = '%-7s %8d %10.6f  %10.8f %2s %14.8f %14.8f\n'
            for row in df.itertuples():
                writer.write(row_fmt % (row.names, row.atnum, row.mass,
                                        row.charge, row.ptype, row.sigma,
                                        row.epsilon))
                writer.write(row_fmt % (f"d{row.names}", row.atnum, row.mass,
                                        row.charge, row.ptype, 0.0, 0.0))
            writer.write('\n')
            writer.write("[ nonbond_params ]\n")
            writer.write("; i     j       func      sigma         epsilon\n")
            row_fmt = '%-7s %-7s %2d %14.8f %14.8f\n'
            for row1 in df.itertuples():
                for row2 in df.itertuples():
                    sigma = (row1.sigma + row2.sigma) / 2
                    epsilon = math.sqrt(row1.epsilon * row2.epsilon)
                    writer.write(row_fmt % (row1.names, row2.names, 1,
                                            sigma, epsilon))
                    writer.write(row_fmt % (f"d{row1.names}", f"d{row2.names}", 1,
                                            sigma, epsilon))
            writer.write('\n')
        return

    def gen_itps(lig_itp, atomtype_itp, decouple_method):
        print(atomtype_itp)
        if decouple_method == 'ci_no' or decouple_method == 'ci_yes':
            for stage in ['coul', 'vdw', 'rest']:
                out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+stage+'.itp'
                atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+stage+'.itp'
                shutil.copyfile(lig_itp, out_name)
                shutil.copyfile(atomtype_itp, atom_out_name)
        elif decouple_method == 'dummy':
            rest_out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+'rest'+'.itp'
            coul_out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+'coul'+'.itp'
            vdw_out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+'vdw'+'.itp'
            rest_atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+'rest'+'.itp'
            coul_atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+'coul'+'.itp'
            vdw_atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+'vdw'+'.itp'
            shutil.copyfile(lig_itp, rest_out_name)
            shutil.copyfile(atomtype_itp, rest_atom_out_name)
            gen_dummy_atomtypes(atomtype_itp, 'coul', coul_atom_out_name)
            gen_dummy_atomtypes(atomtype_itp, 'vdw', vdw_atom_out_name)
            gen_dummy_itps(lig_itp, 'coul', coul_out_name)
            gen_dummy_itps(lig_itp, 'vdw',  vdw_out_name )
        elif decouple_method == 'dummy-pertuball':
            rest_out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+'rest'+'.itp'
            coul_out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+'coul'+'.itp'
            vdw_out_name = args.lig_itp[:-len(args.lig_itp.split('/')[-1])]+ args.lig_itp.split('/')[-1][:-4]+'_'+'vdw'+'.itp'
            rest_atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+'rest'+'.itp'
            coul_atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+'coul'+'.itp'
            vdw_atom_out_name = args.atomtype_itp[:-len(args.atomtype_itp.split('/')[-1])]+ args.atomtype_itp.split('/')[-1][:-4]+'_'+'vdw'+'.itp'
            gen_dummy_atomtypes(atomtype_itp, 'coul', rest_atom_out_name)
            gen_dummy_atomtypes(atomtype_itp, 'coul', coul_atom_out_name)
            gen_dummy_atomtypes(atomtype_itp, 'coul', vdw_atom_out_name)
            gen_dummy_itps(lig_itp, 'coul', rest_out_name)
            gen_dummy_itps(lig_itp, 'coul', coul_out_name)
            gen_dummy_itps(lig_itp, 'coul', vdw_out_name)

    acceptable_alch_ion_methods = ('co-annihilation', 'charge-transfer')
    acceptable_decouple_args = ('ci_no', 'ci_yes', 'dummy', 'dummy-pertuball')

    args = parse_args()
    check_args(args)

    charges = create_charges()
    alch_ion, ion, ion_mass = get_counter_ion_info(args.lig_charge, args.alch_ion_method)
    
    assert get_system_charge(args.gro, alch_ion, ion, ion_mass) == 0, (f'\n{args.gro} is not neutral.\nThis should not happen' 
        '- check your inputs are correct!\nDid you provide the correct "--lig_charge" value?\n(If they are - my bad.)')
    print('Input system is neutral - nice stuff!')

    for stage in ['coul', 'vdw', 'rest']:
        charge_a, charge_b = add_counter_ion_to_gro(args.gro, stage, args.lig_charge, args.alch_ion_method, alch_ion, ion, ion_mass)
        add_counter_ion_to_top(alch_ion, ion, ion_mass, charge_a,charge_b, stage, f'{args.gro[:-4]}_{stage}.gro', args.top)

    if args.gen_itps == True:
        print(f'Generating ITPs using {args.decouple_method}')
        gen_itps(args.lig_itp, args.atomtype_itp, args.decouple_method)

    # add sorting out mdp files
    # even later add functionality to sort out itp files
    # even laterer add functionality to do analysis for a series
