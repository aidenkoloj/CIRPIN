#!/bin/bash

# Counter for processed files
success_count=0
error_count=0

# Process each FASTA file in current directory
for fasta_file in *.fasta; do
    # Check if any FASTA files exist
    if [ ! -e "$fasta_file" ]; then
        echo "No FASTA files found in current directory"
        exit 1
    fi
    
    # Get base filename without extension
    base_name=$(basename "$fasta_file" .fasta)
    output_file="${base_name}_dup.fasta"
    
    echo "Processing: $fasta_file"
    
    # Read the FASTA file and duplicate the sequence
    header=""
    sequence=""
    
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            # This is a header line
            header=$line
        else
            # This is a sequence line
            sequence+=$line
        fi
    done < "$fasta_file"
    
    # Check if we got both header and sequence
    if [ -z "$header" ] || [ -z "$sequence" ]; then
        echo "  ERROR: Invalid FASTA format or empty file"
        ((error_count++))
        continue
    fi
    
    # Modify header to add _dup
    new_header=$(echo "$header" | sed 's/^>\(.*\)$/>\1_dup/')
    
    # Duplicate the sequence
    duplicated_sequence="${sequence}${sequence}"
    
    # Write to output file
    echo "$new_header" > "$output_file"
    echo "$duplicated_sequence" >> "$output_file"
    
    echo "  SUCCESS: Created $output_file"
    ((success_count++))
done

# Print summary
echo ""
echo "=== Duplication Summary ==="
echo "Successful: $success_count"
echo "Failed: $error_count"