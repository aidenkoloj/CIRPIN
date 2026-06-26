#!/bin/bash
BASE_DIR="${1:-.}"
OUTPUT_DIR="/home/gridsan/akolodziej/ECOD_2/ECOD_results/groups_with_cirpin_hits"
N_DIRS=1000
mkdir -p "$OUTPUT_DIR"
SUMMARY_FILE="${OUTPUT_DIR}/scan_summary.txt"
HITS_FILE="${OUTPUT_DIR}/groups_with_hits.txt"
echo "Scanning first $N_DIRS subdirs (by time) in: $BASE_DIR"
echo "---"
copied=0
skipped_empty=0
mapfile -t subdirs < <(
    find "$BASE_DIR" -maxdepth 1 -mindepth 1 -type d \
        ! -path "$OUTPUT_DIR" \
        -printf '%T@ %p\n' \
    | sort -n \
    | head -n "$N_DIRS" \
    | awk '{print $2}'
)
echo "Found ${#subdirs[@]} subdirectories to check"
echo "---"
# Write headers
echo "# All groups checked (${#subdirs[@]} total) — $(date)" > "$SUMMARY_FILE"
printf "%-20s %-12s %s\n" "# group_name" "ecod_group" "hierarchy" >> "$SUMMARY_FILE"
echo "# Groups with CIRPIN hits — $(date)" > "$HITS_FILE"
printf "%-20s %-12s %-8s %s\n" "# group_name" "ecod_group" "n_pairs" "hierarchy" >> "$HITS_FILE"

get_ecod_group() {
    local group_num="$1"
    local group_file="/home/gridsan/akolodziej/ECOD_2/h_groups/group_${group_num}.txt"
    if [[ -f "$group_file" ]]; then
        awk -F'\t' 'NR==1 {print $4}' "$group_file" | cut -d'.' -f1-2
    else
        echo "N/A"
    fi
}

get_ecod_hierarchy() {
    local group_num="$1"
    local group_file="/home/gridsan/akolodziej/ECOD_2/h_groups/group_${group_num}.txt"
    if [[ -f "$group_file" ]]; then
        # Fields 9, 10, 11 (1-indexed in awk) = architecture, x-group, h-group
        awk -F'\t' 'NR==1 {print $9 "\t" $10 "\t" $11}' "$group_file"
    else
        echo "N/A"
    fi
}

for subdir in "${subdirs[@]}"; do
    group_name=$(basename "$subdir")
    group_num="${group_name#group_}"
    filter_file="${subdir}/${group_name}_filtered_results.txt"

    ecod_group=$(get_ecod_group "$group_num")
    hierarchy=$(get_ecod_hierarchy "$group_num")

    # Log every group we attempt
    printf "%-20s %-12s %s\n" "$group_name" "$ecod_group" "$hierarchy" >> "$SUMMARY_FILE"

    [[ -f "$filter_file" ]] || continue

    # Count non-header, non-blank lines
    n_pairs=$(grep -v '^struct1' "$filter_file" | grep -c '\S' || true)

    if [[ "$n_pairs" -gt 0 ]]; then
        cp "$filter_file" "$OUTPUT_DIR/"
        printf "%-20s %-12s %-8s %s\n" "$group_name" "$ecod_group" "${n_pairs}p" "$hierarchy" >> "$HITS_FILE"
        echo "  COPIED: $group_name  (ecod: $ecod_group, $n_pairs pairs) | $hierarchy"
        ((copied++))
    else
        ((skipped_empty++))
    fi
done

echo "" >> "$SUMMARY_FILE"
echo "# Checked: ${#subdirs[@]}  |  Hits: $copied  |  Empty: $skipped_empty" >> "$SUMMARY_FILE"
echo ""
echo "=== Done ==="
echo "  Copied:        $copied"
echo "  Header-only:   $skipped_empty"
echo ""
echo "  Summary:       $SUMMARY_FILE"
echo "  Hits list:     $HITS_FILE"