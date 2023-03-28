import re
from typing import List


class ParsedXVG:
    title: str
    x_label: str
    y_label: str
    col_names: List[str]
    data_rows: List[dict]


def filter_lines(line: str):
    if len(line) == 0:
        return False
    if line[0].strip() == "#":
        return False
    return True


def xvg_to_struct(filename: str) -> ParsedXVG:
    with open(filename) as inp:
        data = list(filter(filter_lines, inp.read().splitlines(keepends=False)))
    res = ParsedXVG()
    legend = []
    for line in filter(lambda x: x[0] == "@", data):
        r = re.match(r'@\s+xaxis\s+label\s+"(?P<lbl>.+)"', line)
        if r:
            res.x_label = r.group('lbl')
            continue
        r = re.match(r'@\s+yaxis\s+label\s+"(?P<lbl>.+)"', line)
        if r:
            res.y_label = r.group('lbl')
            continue
        r = re.match(r'@\s+title\s+"(?P<lbl>.+)"', line)
        if r:
            res.title = r.group('lbl')
            continue
        r = re.match(r'@\s+s(?P<no>\d+)\slegend\s"(?P<lbl>.+)"', line)
        if r:
            legend.append({
                "no": int(r.group('no')),
                "name": r.group('lbl')
            })
        if len(legend) < 1:
            legend.append({
                'no': 1,
                'name': res.y_label
            })
    legend.sort(key=lambda x: x["no"])
    res.col_names = [res.x_label] + [x["name"] for x in legend]
    res.data_rows = []
    for line in filter(lambda x: x[0] != "@", data):
        split = list(filter(lambda x: x != '', line.split()))
        row = {}
        for idx, val in enumerate(split):
            row[res.col_names[idx]] = val
        res.data_rows.append(row)
    return res


if __name__ == "__main__":
    GRAPH = "/media/gluck/Fastdata/big_md/2hzy_charmm36/rmsd.xvg"

    print(max([x['RMSD (nm)'] for x in xvg_to_struct(GRAPH).data_rows]))