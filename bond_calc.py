from pymol import cmd
from typing import List
from dataclasses import dataclass
import json


@dataclass
class Bond:
    name: str
    start_atom: str
    end_atom: str
    start_amino_type: str
    end_amino_type: str
    use_ca: bool
    single_atom: bool
    dist_min: float
    dist_max: float
    inner_color: str
    outer_color: str
    start_amino_list: List[str]
    end_amino_list: List[str]


def get_sec_structs(values):
    print(values)
    for line in values:
        if line.get('sec_letter', 'T') in ('T', 'B', 'G'):
            line['sec_letter'] = 'C'
    structs = []
    curr_struct = None
    for idx, line in enumerate(values):
        if idx == 0:
            curr_struct = {
                "sec_letter": line["sec_letter"],
                "start": int(line["resi"])
            }
        else:
            if line["sec_letter"] != curr_struct["sec_letter"]:
                curr_struct["end"] = int(line["resi"]) - 1
                structs.append(curr_struct)
                curr_struct = {
                    "sec_letter": line["sec_letter"],
                    "start": int(line["resi"])
                }
    curr_struct["end"] = int(values[-1]["resi"])
    structs.append(curr_struct)
    return structs


bonds_data = '[{"name": "disulph", "start_atom": "S", "end_atom": "S", "start_amino_type": "cys", "end_amino_type": "cys", "single_atom": true, "use_ca": false, "dist_min": 1.8, "dist_max": 3.0, "inner_color": "0xAAAAAA", "outer_color": "grey", "start_amino_list": ["CYS"], "end_amino_list": ["CYS"]}, {"name": "hydrogen", "start_atom": "H", "end_atom": "O", "start_amino_type": "d", "end_amino_type": "a", "single_atom": true, "use_ca": false, "dist_min": 1.5, "dist_max": 3.38, "inner_color": "red", "outer_color": "blue", "start_amino_list": ["ARG", "ASN", "GLN", "HIS", "LYS", "SER", "THR", "TRP", "TYR"], "end_amino_list": ["ASN", "ASP", "GLN", "GLU", "HIS", "SER", "THR", "TYR"]}, {"name": "hydrogen", "start_atom": "O", "end_atom": "H", "start_amino_type": "a", "end_amino_type": "d", "single_atom": true, "use_ca": false, "dist_min": 1.5, "dist_max": 3.38, "inner_color": "lightblue", "outer_color": "blue", "start_amino_list": ["ASN", "ASP", "GLN", "GLU", "HIS", "SER", "THR", "TYR"], "end_amino_list": ["ARG", "ASN", "GLN", "HIS", "LYS", "SER", "THR", "TRP", "TYR"]}, {"name": "hydrogen", "start_atom": "H", "end_atom": "N", "start_amino_type": "d", "end_amino_type": "a", "single_atom": true, "use_ca": false, "dist_min": 1.5, "dist_max": 3.38, "inner_color": "lightblue", "outer_color": "blue", "start_amino_list": ["ARG", "ASN", "GLN", "HIS", "LYS", "SER", "THR", "TRP", "TYR"], "end_amino_list": ["ASN", "ASP", "GLN", "GLU", "HIS", "SER", "THR", "TYR"]}, {"name": "hydrogen", "start_atom": "N", "end_atom": "H", "start_amino_type": "a", "end_amino_type": "d", "single_atom": true, "use_ca": false, "dist_min": 1.5, "dist_max": 3.38, "inner_color": "lightblue", "outer_color": "blue", "start_amino_list": ["ASN", "ASP", "GLN", "GLU", "HIS", "SER", "THR", "TYR"], "end_amino_list": ["ARG", "ASN", "GLN", "HIS", "LYS", "SER", "THR", "TRP", "TYR"]}, {"name": "hydrophobic", "start_atom": "", "end_atom": "", "start_amino_type": "hp", "end_amino_type": "hp", "single_atom": false, "use_ca": true, "dist_min": 3.8, "dist_max": 9.5, "inner_color": "red", "outer_color": "red", "start_amino_list": ["ALA", "CYS", "ILE", "LEU", "MET", "PHE", "TRP", "VAL"], "end_amino_list": ["ALA", "CYS", "ILE", "LEU", "MET", "PHE", "TRP", "VAL"]}]'
bond_names = ['disulph', 'hydrogen', 'hydrophobic']

bond_types = [Bond(**x) for x in json.loads(bonds_data)]


def find_bonds(
        motif_sele: str,
        ssec: str = '',
        need_outer=True,
        draw=True,
        remove_surround=False,
        quiet=False):
    bonds = {}
    cmd.select('motif', motif_sele)
    for n in bond_names:
        bonds[f'inner_{n}'] = 0
        if need_outer:
            bonds[f'outer_chain_{n}'] = 0
            bonds[f'outer_interchain_{n}'] = 0
            bonds[f'outer_{n}'] = 0
    motif_atoms = []
    cmd.iterate('%motif and n. CA', lambda atom: motif_atoms.append({'chain': atom.chain, 'resi': int(atom.resi)}))
    ch_list = list(set([x['chain'] for x in motif_atoms]))
    start = motif_atoms[0]['resi']
    if len(ch_list) != 1:
        raise Exception("Bad selection! Multiple chains or no chain!")
    chain = ch_list[0]
    if ssec != '':
        sstruct = []
        cmd.iterate(f'chain {chain} and n. CA', lambda atom: sstruct.append({'resi': int(atom.resi)}))
        for idx, letter in enumerate(list(ssec)):
            sstruct[idx]['sec_letter'] = letter
        sstructs = get_sec_structs(sstruct)
    else:
        sstructs = None
    if remove_surround:
        surround_size = max([x.dist_max for x in bond_types]) * 2
        cmd.remove(f'not (all within {surround_size} of %motif)')
    for bond in bond_types:
        if bond.use_ca:
            start_atom_sele = "n. CA"
            end_atom_sele = "n. CA"
        else:
            start_atom_sele = f"elem {bond.start_atom}"
            end_atom_sele = f"elem {bond.end_atom}"
        start_sele = f"motif and {start_atom_sele} and resn {'+'.join(bond.start_amino_list)}"
        end_atom_sele = f"{end_atom_sele} and resn {'+'.join(bond.end_amino_list)}"
        already_bonded = {k: [] for k in bond_names}

        def iterate_atom(atom):
            if sstructs is not None:
                try:
                    sec_struct = filter(
                        lambda x: x['start'] <= atom.index <= x['end'],
                        sstructs).__next__()
                    if abs(int(atom.resi) - sec_struct['start']) < 3:
                        skip_start = int(atom.resi) - 2
                    else:
                        skip_start = sec_struct['start']
                    if abs(sec_struct['end'] - int(atom.resi)) < 3:
                        skip_end = int(atom.resi) + 2
                    else:
                        skip_end = sec_struct['end']
                except:
                    skip_start = int(atom.resi) - 2
                    skip_end = int(atom.resi) + 2
            else:
                skip_start = int(atom.resi) - 2
                skip_end = int(atom.resi) + 2
            skip_sele = f"chain {chain} and resi {skip_start}-{skip_end}"
            inner_skip = f"resi {start}-{atom.resi}"
            inner_end_sele = f"(%motif and not {inner_skip}) and {end_atom_sele} and not ({skip_sele})"
            outer_chain_end_sele = f"(not motif) and ({end_atom_sele}) and (chain {chain}) and not ({skip_sele})"
            outer_inter_chain_end_sele = f"{end_atom_sele} and not chain {chain}"
            center_sele = f"index {atom.index}"
            cmd.select("distance",
                       f"(({end_atom_sele}) within {bond.dist_max} of {center_sele}) and not"
                       f" (({end_atom_sele}) within {bond.dist_min} of {center_sele})"
                       )
            distances = f"%distance and ({{0}})"
            inner_sele = distances.format(inner_end_sele)

            def iter_t2(atm, bond_type_name, is_inner=True):
                if bond.single_atom:
                    atom_idf = f"{atom.resi}_{atom.chain}"
                    atm_idf = f"{atm.resi}_{atm.chain}"
                    if atom_idf in already_bonded[bond.name] or atm_idf in already_bonded[bond.name]:
                        return
                    already_bonded[bond.name].append(atom_idf)
                    already_bonded[bond.name].append(atm_idf)
                bonds[bond_type_name] += 1
                if draw:
                    name = f"{bond.name}_{atom.resi}_{atm.resi}"
                    cmd.distance(name, f"idx {atom.index}", f"idx {atm.index}")
                    cmd.color(bond.inner_color if is_inner else bond.outer_color, name)

            cmd.iterate(inner_sele, lambda atm: iter_t2(atm, f'inner_{bond.name}'))
            if need_outer:
                cmd.iterate(
                    distances.format(outer_chain_end_sele),
                    lambda atm: iter_t2(atm,
                                        f'outer_chain_{bond.name}',
                                        False)
                )
                cmd.iterate(
                    distances.format(outer_inter_chain_end_sele),
                    lambda atm: iter_t2(atm,
                                        f'outer_interchain_{bond.name}',
                                        False)
                )

        cmd.iterate(start_sele, lambda x: iterate_atom(x))
    if draw:
        cmd.hide('labels')
    if need_outer:
        for n in bond_names:
            bonds[f"outer_{n}"] = bonds[f"outer_chain_{n}"] + bonds[f"outer_interchain_{n}"]
            try:
                bonds[f"{n}_ratio"] = bonds[f"outer_{n}"] / bonds[f"inner_{n}"]
            except ZeroDivisionError:
                bonds[f"{n}_ratio"] = None
            try:
                bonds[f"{n}_chain_ratio"] = bonds[f"outer_interchain_{n}"] / bonds[f"outer_chain_{n}"]
            except ZeroDivisionError:
                bonds[f"{n}_chain_ratio"] = None
    if not quiet:
        print(json.dumps(bonds, indent=2))
    return bonds


cmd.extend('find_bonds', find_bonds)