print("starting...")
import csv
import logging
import datetime
import os

from wrappers import get_all_by_apfid
from apfid import Apfid

SET_NAME = os.environ.get('SETNAME', 'sh3_set_6nmwA')
print(f"opening set {SET_NAME}")
with open(f"/data/sets/{SET_NAME}.csv") as in_file, open(f"/data/res/{SET_NAME}_results_{datetime.datetime.now()}.csv", 'w') as out_file:
    reader = csv.DictReader(in_file)
    writer = csv.DictWriter(out_file, (
        'apfid',
        'desc',
        'src',
        'inner_hydrophobic',
        'outer_hydrophobic',
        'inner_hydrogen',
        'outer_hydrogen',
        'hydrophobic_ratio',
        'hydrogen_ratio',
        'inner_disulph',
        'outer_disulph',
        'disulph_ratio',
        'hydrophobic_chain_ratio',
        'hydrogen_chain_ratio',
        'disulph_chain_ratio',
        'outer_chain_hydrophobic',
        'outer_interchain_hydrophobic',
        'outer_chain_disulph',
        'outer_interchain_disulph',
        'outer_chain_hydrogen',
        'outer_interchain_hydrogen',
        'total_len',
        'chain_len',
        'motif_len',
        'chains',
        'mot_total_ratio',
        'mot_chain_ratio',
        'error',
        'sasa_motif', 'rg', 'sasa_int', 'fasta', 'stride'
    ))
    writer.writeheader()
    out_file.flush()
    for line in reader:
        if line['res'] != 'true':
            continue
        apfid = Apfid(
            line['pdb_id'],
            line['chain'],
            line['start'],
            line['end']
        )
        print(f'writing {apfid.upper()}')
        row = {
            "apfid": apfid.upper(),
            "src": line['src'],
            "desc": line['desc']
        }
        try:
            row.update(get_all_by_apfid(
                apfid,
                '/pdbs',
                use_sec=False
            ))
        except Exception as e:
            logging.exception(e)
            logging.error(line['pdb_id'])
            row['error'] = str(e)
        writer.writerow(row)
        out_file.flush()