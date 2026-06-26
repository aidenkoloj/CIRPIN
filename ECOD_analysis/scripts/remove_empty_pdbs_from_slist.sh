#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --mem=16G

LOG_FILE="structure_list_cleanup.log"
> "$LOG_FILE"  # Clear the log file

echo "Starting structure list validation..." | tee -a "$LOG_FILE"
echo "$(date)" >> "$LOG_FILE"
echo "---" >> "$LOG_FILE"

for structure_list in /home/gridsan/akolodziej/ECOD_2/structure_lists/structure_lists_greater_than_10/structure_list_*.txt; do
    echo "Checking $structure_list..."
    
    needs_resave=false
    temp_file=$(mktemp)
    removed_count=0
    
    while IFS= read -r line; do
        # Extract the path (first column)
        pdb_path=$(echo "$line" | awk '{print $1}')
        
        # Check if file exists and is not empty (using -s flag)
        if [ -s "$pdb_path" ]; then
            # File exists and is not empty, keep it
            echo "$line" >> "$temp_file"
        else
            # File is empty or doesn't exist
            needs_resave=true
            removed_count=$((removed_count + 1))
            echo "  Removed: $pdb_path" | tee -a "$LOG_FILE"
        fi
    done < "$structure_list"
    
    # If modifications were needed, resave the file
    if [ "$needs_resave" = true ]; then
        mv "$temp_file" "$structure_list"
        echo "✓ $structure_list RESAVED - Removed $removed_count empty entries" | tee -a "$LOG_FILE"
    else
        rm "$temp_file"
        echo "✓ $structure_list - All entries valid" | tee -a "$LOG_FILE"
    fi
    echo "---" >> "$LOG_FILE"
done

echo ""
echo "Cleanup complete. Summary saved to: $LOG_FILE"
