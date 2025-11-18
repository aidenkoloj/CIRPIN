### Generate a list of all TED domains that have a PDZ CATH label (2.30.42.10)
### from the AFDB representative clusters.

cluster_fp = '/home/gridsan/akolodziej/TED/ted_365_chunks/cluster_reps_ted_365m.domain_summary.cath.globularity.taxid.tsv'

TED_pdz_domains = []

# CATH ID to search for
CATH_ID = "2.30.42.10"

# List to store matching AF TED IDs
af_ted_ids = []

# Read the file and extract matching IDs
with open(cluster_fp, 'r') as f:
    for line in f:
        if CATH_ID in line:
            # Split by tab and get the first field (AF TED ID)
            af_ted_id = line.split('\t')[0]
            af_ted_ids.append(af_ted_id)

# Print results
print(f"Found {len(af_ted_ids)} entries with CATH ID {CATH_ID}\n")

with open('TED_PDZ_domains.txt', 'w') as f:
    f.write('\n'.join(af_ted_ids))
