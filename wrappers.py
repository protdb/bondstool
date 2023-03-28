import os
import gzip
from pymol import cmd
from bond_calc import find_bonds, bond_names
from stride import get_parsed_stride
from Bio.PDB import PDBParser, PDBIO
from rgyrate import rgyrate
import requests
from apfid import Apfid
from typing import Union

PDB_URL = "https://files.rcsb.org/download/"


def unzip_request(url):
    return gzip.decompress(requests.get(url).content)


parser = PDBParser(QUIET=True)


def get_geometry(file_path, apfid: Union[Apfid, None] = None, sstruc='', need_outer=True, replace_h=False):
    cmd.reinitialize()
    cmd.load(file_path)
    cmd.remove('not polymer')
    if replace_h:
        cmd.remove('hydrogens')
        cmd.h_add()

    if apfid is not None:
        if apfid.start != apfid.end:
            motif_sele = f"resi {apfid.start}-{apfid.end} and chain {apfid.chain_id}"
        else:
            motif_sele = f"chain {apfid.chain_id}"
        res = {
            "total_len": int(cmd.count_atoms("n. CA")),
            "chain_len": int(cmd.count_atoms(f"chain {apfid.chain_id} and n. CA")),
            "chains": len(cmd.get_chains()),
        }
    else:
        motif_sele = 'all'
        res = {}
    res.update({
        "motif_len": int(cmd.count_atoms(f"({motif_sele}) and n. CA")),
        "sasa_int": float(cmd.get_area(motif_sele)),
        "rg": float(rgyrate(motif_sele))
    })
    res.update(find_bonds(
        motif_sele,
        sstruc,
        need_outer,
        False, True, True
    ))
    cmd.remove('not %motif')
    res['sasa_motif'] = cmd.get_area()
    if need_outer:
        for n in bond_names:
            res[f"outer_{n}"] = res[f"outer_chain_{n}"] + res[f"outer_interchain_{n}"]
            try:
                res[f"{n}_ratio"] = res[f"outer_{n}"] / res[f"inner_{n}"]
            except ZeroDivisionError:
                res[f"{n}_ratio"] = None
            try:
                res[f"{n}_chain_ratio"] = res[f"outer_interchain_{n}"] / res[f"outer_chain_{n}"]
            except ZeroDivisionError:
                res[f"{n}_chain_ratio"] = None
    return res


def get_all_by_apfid(apfid, pdb_path='pdbs', use_sec=True, need_outer=True):
    apfid = Apfid(apfid=apfid)
    if not os.path.exists(pdb_path):
        os.mkdir(pdb_path)
    target_file_path = os.path.join(pdb_path, apfid.experiment_id.lower()+'.pdb')
    model0_path = os.path.join(pdb_path, apfid.experiment_id.lower()+'_0.pdb')
    if not os.path.exists(target_file_path):
        with open(target_file_path, 'w') as target_file:
            target_file.write(unzip_request(f"{PDB_URL}{apfid.experiment_id.upper()}.pdb.gz").decode('utf8'))
    if not os.path.exists(model0_path):
        model = parser.get_structure(apfid.experiment_id, target_file_path).get_models().__next__()
        pdb_io = PDBIO()
        pdb_io.save(model0_path, model)
    if use_sec:
        sstruc = ''.join([x['sec_letter'] for x in get_parsed_stride(model0_path) if x['chain'] == apfid.chain_id])
    else:
        sstruc = ''
    return get_geometry(
        model0_path,
        apfid,
        sstruc,
        need_outer,
        replace_h=False
    )


def get_from_md(src_dir, top_file='md.gro', traj_file='md.xtc', src_file='conf.pdb', fragment='', skip=100):
    if fragment != '':
        apfid=Apfid(apfid=fragment)
        need_outer = True
    else:
        chain_id = parser.get_structure('UNK', os.path.join(src_dir, src_file)).get_chains().__next__().id
        apfid = Apfid('UNKN', chain_id)
        need_outer = False
    cmd.reinitialize()
    cmd.load(os.path.join(src_dir, top_file), quiet=True)
    cmd.load_traj(os.path.join(src_dir, traj_file))
    cmd.remove('not polymer')
    cmd.alter('(all)', f'chain="{apfid.chain_id}"')
    states = cmd.count_states()
    tmpdir = os.path.join(src_dir, 'tmp_traj')
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    for state in range(1, states, skip):
        cmd.save(os.path.join(tmpdir, f"{state:06d}.pdb"), state=state)
    intact_params = get_geometry(
        os.path.join(src_dir, src_file),
        apfid,
        need_outer=need_outer,
        replace_h=False
        )
    cmd.reinitialize()
    files = os.listdir(tmpdir)
    files.sort()
    results = []
    for file in files:
        fname = os.path.join(tmpdir, file)
        print(fname)
        bonds = {
            'state': file.split('.')[0]
        }
        bonds.update(get_geometry(fname, apfid, replace_h=False, need_outer=need_outer))
        results.append(bonds)
        os.remove(fname)
    os.rmdir(tmpdir)
    return intact_params, results
