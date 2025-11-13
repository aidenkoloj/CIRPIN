import pickle
import sys

def parse_scope_file(filename):
    """
    Parse SCOPe des file and create a dictionary mapping domain names to SCCSs.
    
    Args:
        filename: Path to the SCOPe des file
        
    Returns:
        Dictionary with domain names as keys and SCCSs as values
    """
    scope_dict = {}
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Skip comment lines
                if line.startswith('#'):
                    continue
                
                # Skip empty lines
                if not line.strip():
                    continue
                
                parts = line.split()
                
                # We need at least 5 fields: id, type, sccs, domain_name, pdb_info
                if len(parts) < 5:
                    continue
                
                # Extract SCCS and domain name
                sccs = parts[2]
                domain_name = parts[3]
                
                # Only include entries with valid domain names (not '-')
                if domain_name != '-' and domain_name.startswith('d'):
                    scope_dict[domain_name] = sccs
        
        return scope_dict
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)


def save_pickle(data, original_filename):
    """
    Save dictionary to pickle file with similar name to original.
    
    Args:
        data: Dictionary to save
        original_filename: Original filename to base pickle name on
    """
    # Create pickle filename by replacing extension
    if '.' in original_filename:
        base_name = original_filename.rsplit('.', 1)[0]
    else:
        base_name = original_filename
    
    pickle_filename = f"{base_name}_sid_sccs_dictionary.pkl"
    
    with open(pickle_filename, 'wb') as f:
        pickle.dump(data, f)
    
    print(f"Dictionary saved to: {pickle_filename}")
    print(f"Total entries: {len(data)}")
    return pickle_filename


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <scope_file>")
        print("Example: python script.py dir.des.scope.2.08-2023-01-06.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Parse the file
    print(f"Parsing {input_file}...")
    scope_dict = parse_scope_file(input_file)
    
    # Save to pickle
    save_pickle(scope_dict, input_file)
    
    # Print some sample entries
    print("\nSample entries:")
    for i, (domain, sccs) in enumerate(list(scope_dict.items())[:5]):
        print(f"  {domain} -> {sccs}")