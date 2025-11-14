## This script filters the TED domain annotation data for only domains occuring in the AFDB clustered representatives (https://www.nature.com/articles/s41586-023-06510-w)


afdb_cluster_reps = '/home/gridsan/akolodziej/Foldcomp/afdb_cluster_reps/afdb_rep_v4.lookup'

print("Loading lookup IDs into set...")
lookup_ids = set()
with open(afdb_cluster_reps) as f:
    for line in f:
        # Adjust based on your lookup file format
        # If ID is first column:
        lookup_id = line.split()[1]
        lookup_ids.add(lookup_id)
        

#*************************************************************************
# File obtained from TED Zenodo; large file, too big to upload to Github #
#**************************************************************************
ted_domains = f"ted_365m.domain_summary.cath.globularity.taxid.tsv"


ted_domains_afdb_cluster_reps = f"CIRPIN/AFDB_ClustR/cluster_reps_ted_365m.domain_summary.cath.globularity.taxid.tsv"


matched_count = 0
total_count = 0

print("Processing TED domains...")
with open(ted_domains) as f_in, open(ted_domains_afdb_cluster_reps, 'w') as f_out:
    for line in f_in:        
        cols = line.rstrip('\n').split('\t')  # Use '\t' for TSV files
        
        ted_id = cols[0]
        afdb_id = ted_id.split("_TED")[0]
        
        if afdb_id in lookup_ids:
            # Write the full line to output file
            f_out.write(line)
            matched_count += 1
        
        total_count += 1
        
        # Progress indicator
        if total_count % 100_000 == 0:
            print(f'Processed {total_count:,} lines, found {matched_count:,} matches...')

print(f"\nComplete!")
print(f"Total lines processed: {total_count:,}")
print(f"Matches found: {matched_count:,}")
print(f"Output written to: {ted_domains_afdb_cluster_reps}")