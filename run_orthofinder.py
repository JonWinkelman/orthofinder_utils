import subprocess

def run_orthofinder(fasta_dir, **kwargs):
    """
    Run OrthoFinder with user-defined options.
    
    Parameters:
        fasta_dir (str): Path to the directory containing FASTA proteomes.
        **kwargs: Additional OrthoFinder options as key-value pairs.

    Example Usage:
        run_orthofinder("./data/proteomes", t=24, s="rooted_tree.nwk", M="msa")
    """
    
    # Base command
    cmd = ["orthofinder", "-f", fasta_dir]
    
    # Process keyword arguments
    for key, value in kwargs.items():
        key = f"-{key}" if not key.startswith("--") else key  # Ensure proper flag format
        
        if isinstance(value, bool):  # Handle boolean flags (e.g., --fewer-files)
            if value:
                cmd.append(key)
        else:
            cmd.extend([key, str(value)])  # Add key-value pairs
    
    # Run the command
    print(f"Running: {' '.join(cmd)}")  # Print for debugging
    subprocess.run(cmd, check=True)

