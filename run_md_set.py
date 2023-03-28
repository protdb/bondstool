import os
import json
import logging
from statistics import mean
from pandas import json_normalize

from apfid import Apfid
from wrappers import get_from_md
from parse_xvg_graph import xvg_to_struct

SKIP = 50
ROOT_DIR = '/data'

for folder in filter(lambda x: os.path.isdir(os.path.join(ROOT_DIR, x)), os.listdir(ROOT_DIR)):
    config = json.loads(os.environ.get("APFIDS", "[]"))
    pdb_id = folder[:4].upper()
    try:
        fragment = filter(lambda x: Apfid(apfid=x).experiment_id.upper() == pdb_id, config).__next__()
    except StopIteration:
        fragment = ''
    work_dir = os.path.join(ROOT_DIR, folder, 'intact_md')
    try:
        intact_bonds, results = get_from_md(
            work_dir,
            traj_file='md_fast.xtc',
            fragment=fragment,
            skip=SKIP
        )
        with open(
            os.path.join(work_dir, 'intact_bonds_v2.json'),
            'w'
        ) as out_file:
            out_file.write(json.dumps(intact_bonds))
        df = json_normalize(results)
        df.to_csv(os.path.join(work_dir, 'traj_lines_v2.csv'))
        with open(os.path.join(
                work_dir,
                "avg_values.json"
        ), 'w') as out_file:
            avgs = {
                "avg_rmsd": mean(
                    [float(x['RMSD (nm)']) for x in xvg_to_struct(os.path.join(work_dir, 'rmsd.xvg')).data_rows]
                ),
                "avg_hg_bonds": mean(df['inner_hydrogen']),
                "avg_hp_bonds": mean(df['inner_hydrophobic']),
                "avg_rg": mean(df['rg']),
                "avg_sasa_int": mean(df['sasa_int']),
                "avg_sasa_motif": mean(df['sasa_motif'])
            }
            try:
                avgs.update({
                    "avg_hg_out": mean(df['outer_hydrogen']),
                    "avg_hp_out": mean(df['outer_hydrophobic'])
                })
            except KeyError:
                pass
            out_file.write(json.dumps(avgs, indent=2))

    except Exception as e:
        logging.exception(e)
        continue