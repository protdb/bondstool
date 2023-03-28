import subprocess
import os


def run_stride(fname: str, force: bool = False):
    inp = fname
    outp = fname+'.stride'
    if not os.path.exists(outp) and not force:
        try:
            with open(os.devnull, 'w') as devnul:
                subprocess.check_call(["./stride", '-f' + outp, inp], stderr=devnul, stdout=devnul)
        except subprocess.CalledProcessError:
            return None
    with open(outp) as outfile:
        return outfile.read()


def get_parsed_stride(fname):
    raw = run_stride(fname)
    if raw is None:
        return None
    try:
        asg = filter(lambda x: x[:3] == "ASG", raw.splitlines(keepends=False))
    except AttributeError:
        return None
    parsed = []
    try:
        for line in asg:
            parsed.append({
                "amino": line[5:8],
                "total_no": int(line[10:15]),
                "file_no": int(line[16:20]),
                "chain": line[9],
                "sec_letter": line[24],
                "sec_name": line[25:39].strip()
            })
    except ValueError:
        return None
    return parsed


def get_sec_structures(parsed):
    for line in parsed:
        if line['sec_letter'] == 'T':
            line['sec_letter'] = 'C'
    structs = []
    curr_struct = None
    for idx, line in enumerate(parsed):
        if idx == 0:
            curr_struct = {
                "type": line["sec_name"],
                "sec_letter": line["sec_letter"],
                "start": int(line["total_no"])
            }
        else:
            if line["sec_letter"] != curr_struct["sec_letter"]:
                curr_struct["end"] = int(line["total_no"]) - 1
                structs.append(curr_struct)
                curr_struct = {
                    "type": line["sec_name"],
                    "sec_letter": line["sec_letter"],
                    "start": int(line["total_no"])
                }
    curr_struct["end"] = int(parsed[-1]["total_no"])
    structs.append(curr_struct)
    return structs


class SecondaryStructureException(Exception):
    pass


class StrideFinder:
    file_name = ""
    parsed = []
    sec_structs = []

    def __init__(self, file_name):
        self.file_name = file_name
        self.parsed = get_parsed_stride(file_name)
        if self.parsed is None:
            raise SecondaryStructureException("Empty secondary structure!")
        self.sec_structs = get_sec_structures(self.parsed)

    def get_struct_by_idx(self, idx: int):
        try:
            return filter(
                lambda x: x['start'] <= idx <= x['end'],
                self.sec_structs).__next__()
        except StopIteration:
            raise SecondaryStructureException(f"cannot find sec. struct. in file {self.file_name} for res {idx}!")