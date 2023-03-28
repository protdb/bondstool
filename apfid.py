import logging
from typing import Union


class Apfid:
    apfid: str
    experiment_id: str
    chain_id: str
    start: Union[int, None]
    end: Union[int, None]
    id_type: str

    def __init__(self, experiment_id='', chain_id='', start=0, end=0, apfid=None):
        logging.debug(f"{experiment_id} {chain_id} {start} {end}")
        if apfid is None:
            self.experiment_id = experiment_id
            self.chain_id = chain_id
            if start == end:
                self.start = None
                self.end = None
            else:
                self.start = start
                self.end = end
            self.apfid = self._make_apfid()
        else:
            self.apfid = apfid
            self._parse_apfid()
        self._set_experiment_type()

    def _set_experiment_type(self):
        if len(self.experiment_id) == 4:
            self.id_type = "PDB"
        else:
            self.id_type = "Alphafold"

    def _make_apfid(self, lower=False):
        exp_id = self.experiment_id.lower() if lower else self.experiment_id.upper()
        if self.start == self.end:
            return f"{exp_id}_{self.chain_id}"
        else:
            return f"{exp_id}_{self.chain_id}{self.start}_{self.chain_id}{self.end}"

    def _parse_apfid(self):
        split = self.apfid.split("_")
        self.experiment_id = split[0]
        if len(split) == 2:
            self.chain_id = split[1]
            self.start = None
            self.end = None
        else:
            self.chain_id = split[1][0]
            self.start = int(split[1][1:])
            self.end = int(split[2][1:])

    def __str__(self):
        return self.apfid

    def upper(self):
        return self._make_apfid(lower=False)

    def lower(self):
        return self._make_apfid(lower=True)
