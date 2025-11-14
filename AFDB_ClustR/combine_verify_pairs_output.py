### Script to combine all the outputs from the
# verify putative pairs script
#### Outputs to combine:
# verified_pairs_AFDB_cluster_reps_*
# unqiue_CPs_AFDB_cluster_reps_*
# other_homologous_pairs_AFDB_cluster_reps_*
# false_positive_pairs_AFDB_cluster_reps_*

import os
import glob

top_dir = '/home/gridsan/akolodziej/TED/ted_365_chunks'

# Define the output file patterns to look for
output_patterns = [
    'verified_pairs_AFDB_cluster_reps_*.txt',
    'unique_CPs_AFDB_cluster_reps_*.txt',
    'other_homologous_pairs_AFDB_cluster_reps_*.txt',
    'false_positive_pairs_AFDB_cluster_reps_*.txt'
]

def is_data_line(line, header):
    """Check if a line is a data line (not header or summary stats)"""
    line = line.strip()
    
    # Skip empty lines
    if not line:
        return False
    
    # Skip header lines
    if line == header.strip():
        return False
    
    # Skip summary statistics lines (start with "Number" or "Percent")
    if line.startswith('Number of') or line.startswith('Percent of'):
        return False
    
    # Check if line starts with AF- (typical protein ID format)
    if line.startswith('AF-'):
        return True
    
    return False

# For each output pattern, find all matching files and combine them
for pattern in output_patterns:
    print(f"\nProcessing pattern: {pattern}")
    
    # Find all files matching this pattern
    file_pattern = os.path.join(top_dir, pattern)
    matching_files = sorted(glob.glob(file_pattern))
    
    if not matching_files:
        print(f"  No files found for pattern: {pattern}")
        continue
    
    print(f"  Found {len(matching_files)} files")
    
    # Determine output filename (remove the wildcard)
    output_filename = pattern.replace('_*.txt', '_combined.txt')
    output_path = os.path.join(top_dir, output_filename)
    
    # Combine all files
    header_written = False
    header = None
    total_lines = 0
    
    with open(output_path, 'w') as outfile:
        for file in matching_files:
            try:
                with open(file, 'r') as infile:
                    lines = infile.readlines()
                    
                    if not lines:
                        continue
                    
                    # Get header from first file
                    if header is None:
                        header = lines[0]
                    
                    # Write header only once (from first file)
                    if not header_written:
                        outfile.write(header)
                        header_written = True
                    
                    # Write only valid data lines (skip header and summary stats)
                    file_data_lines = 0
                    for line in lines[1:]:
                        if is_data_line(line, header):
                            outfile.write(line)
                            file_data_lines += 1
                    
                    total_lines += file_data_lines
                    print(f"    Read {file_data_lines} rows from {os.path.basename(file)}")
                    
            except Exception as e:
                print(f"    Error reading {file}: {e}")
    
    print(f"  Saved {total_lines} total rows to {output_filename}")

print("\nAll files combined successfully!")