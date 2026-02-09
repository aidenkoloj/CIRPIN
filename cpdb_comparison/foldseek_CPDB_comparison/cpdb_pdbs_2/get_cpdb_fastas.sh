#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p ./fasta_seqs

# Initialize error log
error_file="./fasta_seqs/conversion_errors.txt"
> "$error_file"

# Counter for successful conversions
success_count=0
error_count=0

# Process each PDB file in current directory
for pdb_file in *.pdb; do
    # Check if any PDB files exist
    if [ ! -e "$pdb_file" ]; then
        echo "No PDB files found in current directory"
        exit 1
    fi
    
    # Get base filename without extension
    base_name=$(basename "$pdb_file" .pdb)
    fasta_file="./fasta_seqs/${base_name}.fasta"
    
    echo "Processing: $pdb_file"
    
    # Run pdb_tofasta command and modify header
    pdb_tofasta "$pdb_file" 2>/dev/null | sed "s/^>.*/>$base_name/" > "$fasta_file"
    
    # Check if conversion was successful
    if [ ! -s "$fasta_file" ]; then
        # File is empty or doesn't exist
        echo "$pdb_file" >> "$error_file"
        echo "  ERROR: Conversion failed (empty output)"
        ((error_count++))
        rm -f "$fasta_file"  # Remove empty file
    else
        echo "  SUCCESS: Created $fasta_file"
        ((success_count++))
    fi
done

# Print summary
echo ""
echo "=== Conversion Summary ==="
echo "Successful: $success_count"
echo "Failed: $error_count"

if [ $error_count -gt 0 ]; then
    echo "Failed files logged in: $error_file"
fi