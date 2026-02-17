#!/usr/bin/env python3
"""
Foldseek Clustering Parameter Sweep

This script runs foldseek easy-cluster with different parameter combinations 
and analyzes how many clusters are produced.
"""

import pandas as pd
import subprocess
import os
import shutil
from pathlib import Path
import matplotlib.pyplot as plt

# ============================================================================
# Setup Parameters
# ============================================================================

# Input directory with PDB files
pdb_dir = '/home/ubuntu/CIRPIN/cpdb_comparison/foldseek_CPDB_comparison/CIRPIN_pdbs/scope40_only'

# Output directory for results
output_base_dir = '/home/ubuntu/CIRPIN/cpdb_comparison/foldseek_CPDB_comparison/CIRPIN_pdbs/scope40_only/parameter_sweep_results'
os.makedirs(output_base_dir, exist_ok=True)

# Temporary directory
tmp_dir = '/home/ubuntu/CIRPIN/cpdb_comparison/foldseek_CPDB_comparison/CIRPIN_pdbs/scope40_only/tmp'
os.makedirs(tmp_dir, exist_ok=True)

# Fixed parameters
e_value = 0.01
coverage = 0.9
cov_mode = 0  # default
tmscore_threshold = 0.5
exact_tmscore = 0  # default (approximate)

# Variable parameters
alignment_types = [0, 1, 2]  # 0: 3di, 1: TM alignment, 2: 3Di+AA
cluster_modes = [0, 1, 2, 3]  # 0: Set-Cover, 1: Connected component, 2,3: CDHIT variants

# Labels for clarity
alignment_labels = {0: '3di', 1: 'TM', 2: '3Di+AA'}
cluster_labels = {0: 'SetCover', 1: 'ConnectedComp', 2: 'CDHIT_v1', 3: 'CDHIT_v2'}

# ============================================================================
# Run Foldseek Clustering
# ============================================================================

results = []

for align_type in alignment_types:
    for clust_mode in cluster_modes:
        # Create descriptive output prefix
        output_prefix = f"align{alignment_labels[align_type]}_clust{cluster_labels[clust_mode]}_e{e_value}_c{coverage}"
        output_path = os.path.join(output_base_dir, output_prefix)
        
        print(f"\n{'='*60}")
        print(f"Running: {output_prefix}")
        print(f"{'='*60}")
        
        # Build foldseek command
        cmd = [
            'foldseek', 'easy-cluster',
            pdb_dir,
            output_path,
            tmp_dir,
            '--alignment-type', str(align_type),
            '--cluster-mode', str(clust_mode),
            '-e', str(e_value),
            '-c', str(coverage),
            '--cov-mode', str(cov_mode),
            '--tmscore-threshold', str(tmscore_threshold),
            '--exact-tmscore', str(exact_tmscore)
        ]
        
        print(f"Command: {' '.join(cmd)}")
        
        # Run foldseek
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print("Success!")
            if result.stdout:
                print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error running foldseek: {e}")
            print(f"stderr: {e.stderr}")
            results.append({
                'alignment_type': align_type,
                'alignment_label': alignment_labels[align_type],
                'cluster_mode': clust_mode,
                'cluster_label': cluster_labels[clust_mode],
                'e_value': e_value,
                'coverage': coverage,
                'cov_mode': cov_mode,
                'tmscore_threshold': tmscore_threshold,
                'num_clusters': 'ERROR',
                'output_prefix': output_prefix
            })
            # Clean up tmp directory after error (foldseek bug workaround)
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
            os.makedirs(tmp_dir, exist_ok=True)
            continue
        
        # Count clusters from the output TSV file
        cluster_file = f"{output_path}_cluster.tsv"
        
        if os.path.exists(cluster_file):
            try:
                # Read cluster file
                fs_cluster = pd.read_csv(cluster_file, sep='\t', header=None)
                fs_cluster.columns = ['Cluster', 'Structure']
                
                # Count unique clusters
                cluster_dict = fs_cluster.groupby('Cluster')['Structure'].apply(list).to_dict()
                num_clusters = len(cluster_dict)
                
                print(f"Number of clusters: {num_clusters}")
                
                results.append({
                    'alignment_type': align_type,
                    'alignment_label': alignment_labels[align_type],
                    'cluster_mode': clust_mode,
                    'cluster_label': cluster_labels[clust_mode],
                    'e_value': e_value,
                    'coverage': coverage,
                    'cov_mode': cov_mode,
                    'tmscore_threshold': tmscore_threshold,
                    'num_clusters': num_clusters,
                    'output_prefix': output_prefix
                })
            except Exception as e:
                print(f"Error reading cluster file: {e}")
                results.append({
                    'alignment_type': align_type,
                    'alignment_label': alignment_labels[align_type],
                    'cluster_mode': clust_mode,
                    'cluster_label': cluster_labels[clust_mode],
                    'e_value': e_value,
                    'coverage': coverage,
                    'cov_mode': cov_mode,
                    'tmscore_threshold': tmscore_threshold,
                    'num_clusters': 'ERROR_READING',
                    'output_prefix': output_prefix
                })
        else:
            print(f"Warning: Cluster file not found: {cluster_file}")
            results.append({
                'alignment_type': align_type,
                'alignment_label': alignment_labels[align_type],
                'cluster_mode': clust_mode,
                'cluster_label': cluster_labels[clust_mode],
                'e_value': e_value,
                'coverage': coverage,
                'cov_mode': cov_mode,
                'tmscore_threshold': tmscore_threshold,
                'num_clusters': 'FILE_NOT_FOUND',
                'output_prefix': output_prefix
            })
        
        # Clean up tmp directory after each run (foldseek bug workaround)
        print("Cleaning up tmp directory...")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir, exist_ok=True)

print(f"\n\n{'='*60}")
print("All runs completed!")
print(f"{'='*60}")

# ============================================================================
# Create Results Summary
# ============================================================================

# Create DataFrame with results
results_df = pd.DataFrame(results)

# Reorder columns for better readability
column_order = [
    'alignment_label',
    'cluster_label',
    'num_clusters',
    'alignment_type',
    'cluster_mode',
    'e_value',
    'coverage',
    'cov_mode',
    'tmscore_threshold',
    'output_prefix'
]
results_df = results_df[column_order]

# Display results
print("\nResults Summary:")
print(results_df.to_string(index=False))

# ============================================================================
# Save Results to CSV
# ============================================================================

# Save to CSV (no extra dependencies needed)
csv_output = os.path.join(output_base_dir, 'clustering_results_summary.csv')
results_df.to_csv(csv_output, index=False)
print(f"\nResults saved to: {csv_output}")

# ============================================================================
# Pivot Table Analysis
# ============================================================================

# Create pivot table showing clusters by alignment type and cluster mode
pivot_df = results_df.pivot_table(
    values='num_clusters',
    index='cluster_label',
    columns='alignment_label',
    aggfunc='first'
)

print("\nPivot Table: Number of Clusters")
print("Rows: Cluster Mode | Columns: Alignment Type")
print(pivot_df)

# ============================================================================
# Visualization
# ============================================================================

# Filter out error rows for plotting
plot_df = results_df[results_df['num_clusters'].apply(lambda x: isinstance(x, (int, float)))].copy()

if len(plot_df) > 0:
    # Create grouped bar plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Prepare data for grouped bars
    alignment_types_plot = plot_df['alignment_label'].unique()
    cluster_modes_plot = plot_df['cluster_label'].unique()
    
    x = range(len(cluster_modes_plot))
    width = 0.25
    
    for i, align_label in enumerate(alignment_types_plot):
        data = plot_df[plot_df['alignment_label'] == align_label]
        data = data.sort_values('cluster_mode')
        offset = width * (i - 1)
        ax.bar([p + offset for p in x], data['num_clusters'], width, label=align_label)
    
    ax.set_xlabel('Cluster Mode')
    ax.set_ylabel('Number of Clusters')
    ax.set_title('Foldseek Clustering Results: Effect of Alignment Type and Cluster Mode')
    ax.set_xticks(x)
    ax.set_xticklabels(cluster_modes_plot)
    ax.legend(title='Alignment Type')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    fig_output = os.path.join(output_base_dir, 'clustering_results_plot.png')
    plt.savefig(fig_output, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {fig_output}")
    
    # Don't show in script mode (only save)
    # plt.show()
else:
    print("No valid data to plot (all runs had errors)")

print("\n" + "="*60)
print("ANALYSIS COMPLETE")
print("="*60)
