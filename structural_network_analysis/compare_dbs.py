from pathlib import Path
import os
import pickle

def tmscore(q,t, cp=False):
    ''' Run TM-align and get back TM-align score '''
    if cp:
        output = os.popen(f'/home/ubuntu/TM_tools/TMalign {q} {t} -cp')
    else:
        output = os.popen(f'/home/ubuntu/TM_tools/TMalign {q} {t}')
    tms = {"tms":[]}
    parse_float = lambda x: float(x.split("=")[1].split()[0])

    for line in output:
        line = line.rstrip()
        if line.startswith("TM-score"): 
            tms["tms"].append(parse_float(line))
    min_tms = min(tms['tms'])
    return min_tms


CIRPIN_pdbs = [
    p for p in Path("/home/ubuntu/CIRPIN_db/ark_CPDB_analysis/CIRPIN_all_structures").iterdir()
    if p.is_file() and p.suffix == ".pdb"
]

CPDB_pdbs = [
    p for p in Path('/home/ubuntu/CIRPIN_db/ark_CPDB_analysis/CPDB_pdb_structures').iterdir()
    if p.is_file() and p.suffix == ".pdb"
]

hits = []

for i,q in enumerate(CPDB_pdbs):
    print(f'Processing {i} from CPDB')
    for j,t in enumerate(CIRPIN_pdbs):
        print(q,t)
        tm = tmscore(q,t, cp=True)
        if tm > 0.5:
            pair = [q,t,tm]
            hits.append(pair)
            print(f'Found hit! {pair}')

with open('overlap_CIRPIN_scope_CPDB.pkl', "wb") as f:
    pickle.dump(hits, f)