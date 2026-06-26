import json

ecod_dict = {}

with open('ecod.v294.1.F40.names.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            try:
                uid = int(parts[0])
                ecod_domain_id = parts[1]
                ecod_dict[ecod_domain_id] = str(uid).zfill(9)
            except ValueError:
                # Skip non-numeric lines (like the header)
                continue

# Save as JSON
with open('ecod_dict.json', 'w') as f:
    json.dump(ecod_dict, f)

print(f"Saved {len(ecod_dict)} entries")