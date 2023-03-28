from dataclasses import dataclass
import csv
import json
from typing import List, Tuple


@dataclass
class Bond:
    name: str
    start_atom: str
    end_atom: str
    start_amino_type: str
    end_amino_type: str
    single_atom: bool
    use_ca: bool
    dist_min: float
    dist_max: float
    inner_color: str
    outer_color: str
    start_amino_list: List[str]
    end_amino_list: List[str]


def get_bonds() -> Tuple[List[Bond], List[str]]:
    with open("residues.json") as res_in_file:
        residues = json.loads(res_in_file.read())
    with open("bond_types.csv") as bonds_in_file:
        reader = csv.DictReader(bonds_in_file)
        out_bonds = []
        for line in reader:
            line['use_ca'] = bool(int(line['use_ca']))
            line['dist_min'] = float(line['dist_min'])
            line['single_atom'] = bool(int(line['single_atom']))
            line['dist_max'] = float(line['dist_max'])
            line['start_amino_list'] = residues[line['start_amino_type']]
            line['end_amino_list'] = residues[line['end_amino_type']]
            out_bonds.append(Bond(**line))
        out_names = list(set([x.name for x in out_bonds]))
        return out_bonds, out_names


if __name__ == "__main__":
    print(json.dumps([x.__dict__ for x in get_bonds()[0]]))