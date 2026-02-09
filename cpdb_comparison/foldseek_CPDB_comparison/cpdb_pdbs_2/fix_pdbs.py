#### Fix all the PDB files that have repeated chains due to being NMR ensembles, etc

from pathlib import Path
import pandas as pd

def check_and_trim_pdb_inplace(pdb_path):
    """
    Detects residue-number restarts within a chain.
    If detected, trims the file at the restart point
    and overwrites the original PDB.
    """
    pdb_path = Path(pdb_path)

    last_resid = {}
    lines_to_write = []
    repeat_found = False

    with pdb_path.open() as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                resid = int(line[22:26])

                if chain in last_resid and resid < last_resid[chain]:
                    repeat_found = True
                    break

                last_resid[chain] = resid

            lines_to_write.append(line)

    if repeat_found:
        pdb_path.write_text("".join(lines_to_write))
        return True
    else:
        return False


# ---- Fix all PDBs in current directory ----

fixed_files = []
all_files = list(Path.cwd().glob("*.pdb"))

for idx, pdb in enumerate(all_files, 1):
    if check_and_trim_pdb_inplace(pdb):
        fixed_files.append(pdb.name)

    # Print progress every 100 files
    if idx % 100 == 0:
        print(f"Checked {idx} / {len(all_files)} files...")

# Save list of fixed files
if fixed_files:
    df = pd.DataFrame({"fixed_pdb": fixed_files})
    df.to_csv("trimmed_pdbs.csv", index=False)
    print(f"\n📝 Wrote trimmed_pdbs.csv ({len(fixed_files)} files trimmed)")
else:
    print("\n✅ No files required trimming")
