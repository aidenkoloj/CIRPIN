# For each homology group file, extract ECOD domain IDs and create a corresponding structure list mapping them to PDB file paths.

import json
import os
import re

# Load the dictionary
with open('ecod_dict.json', 'r') as f:
    ecod_dict = json.load(f)

h_groups_dir = '/home/gridsan/akolodziej/ECOD_2/h_groups'
structures_dir = '/home/gridsan/groups/solab/ECOD_ark/ecod_structures'

# Loop through all group files
for filename in sorted(os.listdir(h_groups_dir)):
    if filename.startswith('group_') and filename.endswith('.txt'):
        filepath = os.path.join(h_groups_dir, filename)
        
        # Extract group number from filename (e.g., 'group_1.txt' -> '1')
        match = re.search(r'group_(\d+)\.txt', filename)
        if not match:
            continue
        
        group_num = match.group(1)
        output_filename = f'structure_list_{group_num}.txt'
        skipped_filename = f'skipped_structures_{group_num}.txt'
        
        structure_list = []
        skipped_ids = []
        
        # Read each line in the group file
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Extract ecod_domain_id (second column)
                parts = line.split('\t')
                if len(parts) >= 2:
                    ecod_id = parts[1]
                    
                    # Look up in dictionary
                    if ecod_id in ecod_dict:
                        pdb_id = ecod_dict[ecod_id]
                        pdb_path = f'{structures_dir}/{pdb_id}.pdb'
                        structure_list.append(f'{pdb_path} {pdb_id}')
                    else:
                        skipped_ids.append(ecod_id)
        
        # Write structure_list for this group
        with open(output_filename, 'w') as f:
            for line in structure_list:
                f.write(line + '\n')
        
        # Write skipped IDs if any
        if skipped_ids:
            with open(skipped_filename, 'w') as f:
                for ecod_id in skipped_ids:
                    f.write(ecod_id + '\n')
            print(f"Created {output_filename} with {len(structure_list)} entries, {len(skipped_ids)} skipped (see {skipped_filename})")
        else:
            print(f"Created {output_filename} with {len(structure_list)} entries")