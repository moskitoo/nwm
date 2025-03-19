# import os
# import gzip
# import numpy as np
# from collections import defaultdict
# from prettytable import PrettyTable
# import resource
# from multiprocessing import Pool, cpu_count

# # Standard amino acids
# STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')

# def calculate_aa_percentage(sequences, aa):
#     """Calculate the percentage of a specific amino acid in a list of sequences."""
#     total_count = sum(seq.count(aa) for seq in sequences)
#     total_length = sum(len(seq) for seq in sequences)
#     return (total_count / total_length) * 100 if total_length > 0 else 0

# def print_memory_usage():
#     """Print current memory usage."""
#     usage = resource.getrusage(resource.RUSAGE_SELF)
#     print(f"Memory usage: {usage.ru_maxrss / 1024:.2f} MB")

# def read_fasta(file_path):
#     """Read a FASTA file (or .gz) and yield sequences."""
#     try:
#         if file_path.endswith('.gz'):
#             with gzip.open(file_path, 'rt') as file:
#                 sequence = ''
#                 for line in file:
#                     line = line.strip()
#                     if line.startswith('>'):
#                         if sequence:
#                             yield sequence
#                         sequence = ''
#                     else:
#                         sequence += line
#                 if sequence:
#                     yield sequence
#         else:
#             with open(file_path, 'r') as file:
#                 sequence = ''
#                 for line in file:
#                     line = line.strip()
#                     if line.startswith('>'):
#                         if sequence:
#                             yield sequence
#                         sequence = ''
#                     else:
#                         sequence += line
#                 if sequence:
#                     yield sequence
#     except Exception as e:
#         print(f"Error reading {file_path}: {e}")

# def process_file(file_path):
#     """Process a single FASTA file and return amino acid counts and total length."""
#     aa_counts = defaultdict(int)
#     total_length = 0
#     sequences = []

#     for seq in read_fasta(file_path):
#         total_length += len(seq)
#         sequences.append(seq)
#         for aa in seq:
#             if aa in STANDARD_AMINO_ACIDS:
#                 aa_counts[aa] += 1

#     return aa_counts, total_length, sequences

# def bootstrap_iteration(args):
#     """Perform a single bootstrap iteration."""
#     sequences, func, aa, i, n_bootstrap = args
#     print(f"Bootstrap iteration {i+1}/{n_bootstrap}")  # Debugging
#     print_memory_usage()  # Debugging
#     # Resample with replacement
#     resampled_sequences = np.random.choice(sequences, size=len(sequences), replace=True)
#     if aa is not None:
#         return func(resampled_sequences, aa)  # For functions that need `aa`
#     else:
#         return func(resampled_sequences)  # For functions that don't need `aa`

# def bootstrap_error(sequences, func, aa=None, n_bootstrap=3):
#     """Estimate error using bootstrapping with parallel processing."""
#     # Prepare arguments for bootstrap_iteration
#     args = [(sequences, func, aa, i, n_bootstrap) for i in range(n_bootstrap)]

#     # Use multiprocessing to parallelize bootstrapping
#     with Pool() as pool:
#         values = pool.map(bootstrap_iteration, args)

#     return np.std(values)  # Standard deviation as error

# def process_group_parallel(directory, group_name):
#     """Process a group of FASTA files in parallel and calculate averages."""
#     aa_counts = defaultdict(int)
#     total_length = 0
#     all_sequences = []

#     print(f"Processing directory: {directory} for group: {group_name}")
#     file_paths = []
#     for root, _, files in os.walk(directory):
#         for file in files:
#             if file.startswith(group_name) and (file.endswith('.fasta.gz') or file.endswith('.fasta')):
#                 file_paths.append(os.path.join(root, file))

#     # Use multiprocessing to process files in parallel
#     with Pool(cpu_count()) as pool:
#         results = pool.map(process_file, file_paths)

#     for result in results:
#         counts, length, sequences = result
#         for aa, count in counts.items():
#             aa_counts[aa] += count
#         total_length += length
#         all_sequences.extend(sequences)

#     if not all_sequences:
#         print(f"No sequences found for {group_name} in {directory}")
#         return {}, {}, 0, 0

#     # Calculate average aa content
#     aa_content = {aa: (count / total_length) * 100 for aa, count in aa_counts.items()}
#     print(f"Calculating amino acid content errors for {group_name}...")
#     aa_content_error = {}
#     for aa in STANDARD_AMINO_ACIDS:
#         aa_content_error[aa] = bootstrap_error(all_sequences, calculate_aa_percentage, aa)

#     # Calculate average protein length
#     sequence_lengths = [len(seq) for seq in all_sequences]
#     avg_length = np.mean(sequence_lengths)
#     print(f"Calculating protein length error for {group_name}...")
#     length_error = bootstrap_error(sequence_lengths, np.mean, None)

#     return aa_content, aa_content_error, avg_length, length_error

# def main():
#     """Main function to process data and generate tables."""
#     base_dir = os.path.abspath(os.path.dirname(__file__))
#     selected_organisms_dir = os.path.join(base_dir, "selected_organisms")
#     random_selected_dir = os.path.join(base_dir, "random_selected")

#     print(f"Base directory: {base_dir}")
#     print(f"Selected organisms directory: {selected_organisms_dir}")
#     print(f"Random selected directory: {random_selected_dir}")

#     # Process selected organisms (a)
#     print("\nProcessing selected organisms (a)...")
#     aa_content_a, aa_content_error_a, avg_length_a, length_error_a = process_group_parallel(selected_organisms_dir, "Selected Organisms")

#     # Process randomly selected kingdoms (c)
#     print("\nProcessing randomly selected kingdoms (c)...")
#     kingdoms = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
#     aa_content_c = {}
#     aa_content_error_c = {}
#     avg_length_c = {}
#     length_error_c = {}

#     for kingdom in kingdoms:
#         print(f"\nProcessing {kingdom} files in directory: {random_selected_dir}")
#         aa_content_c[kingdom], aa_content_error_c[kingdom], avg_length_c[kingdom], length_error_c[kingdom] = process_group_parallel(random_selected_dir, kingdom)

#     # Create tables using PrettyTable
#     table1 = PrettyTable()
#     table1.field_names = ["Amino Acid", "Average Content (%)", "Error"]
#     for aa in STANDARD_AMINO_ACIDS:
#         table1.add_row([aa, f"{aa_content_a.get(aa, 0):.2f}", f"{aa_content_error_a.get(aa, 0):.2f}"])

#     table2 = PrettyTable()
#     table2.field_names = ["Group", "Average Length", "Error"]
#     table2.add_row(["Selected Organisms", f"{avg_length_a:.2f}", f"{length_error_a:.2f}"])

#     table3 = PrettyTable()
#     table3.field_names = ["Kingdom", "Amino Acid", "Average Content (%)", "Error"]
#     for kingdom in kingdoms:
#         for aa in STANDARD_AMINO_ACIDS:
#             table3.add_row([kingdom, aa, f"{aa_content_c[kingdom].get(aa, 0):.2f}", f"{aa_content_error_c[kingdom].get(aa, 0):.2f}"])

#     table4 = PrettyTable()
#     table4.field_names = ["Kingdom", "Average Length", "Error"]
#     for kingdom in kingdoms:
#         table4.add_row([kingdom, f"{avg_length_c[kingdom]:.2f}", f"{length_error_c[kingdom]:.2f}"])

#     print("\nTable 1: Average Amino Acid Content for Selected Organisms (a)")
#     print(table1)

#     print("\nTable 2: Average Protein Length for Selected Organisms (a)")
#     print(table2)

#     print("\nTable 3: Average Amino Acid Content for Randomly Selected Kingdoms (c)")
#     print(table3)

#     print("\nTable 4: Average Protein Length for Randomly Selected Kingdoms (c)")
#     print(table4)

# if __name__ == '__main__':
#     main()

import os
import gzip
import numpy as np
from collections import defaultdict
from prettytable import PrettyTable
import resource
from multiprocessing import Pool, cpu_count
import random
import gc

# Standard amino acids
STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')

def calculate_aa_percentage(sequences, aa):
    """Calculate the percentage of a specific amino acid in a list of sequences."""
    total_count = sum(seq.count(aa) for seq in sequences)
    total_length = sum(len(seq) for seq in sequences)
    return (total_count / total_length) * 100 if total_length > 0 else 0

def print_memory_usage():
    """Print current memory usage."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    print(f"Memory usage: {usage.ru_maxrss / 1024:.2f} MB")

def read_fasta(file_path, max_sequences=1000):
    """Read a FASTA file (or .gz) and yield sequences with a limit."""
    count = 0
    try:
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = 'rt' if file_path.endswith('.gz') else 'r'
        
        with open_func(file_path, mode) as file:
            sequence = ''
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if sequence:
                        yield sequence
                        count += 1
                        if count >= max_sequences:
                            return
                    sequence = ''
                else:
                    sequence += line
            if sequence:
                yield sequence
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

def process_file(args):
    """Process a single FASTA file and return amino acid counts and total length."""
    file_path, max_sequences, sample_size = args
    aa_counts = defaultdict(int)
    total_length = 0
    sequence_lengths = []
    
    # Sample a limited number of sequences for bootstrapping
    sample_sequences = []
    
    for seq in read_fasta(file_path, max_sequences):
        total_length += len(seq)
        sequence_lengths.append(len(seq))
        
        # Randomly sample for bootstrapping to avoid memory issues
        if random.random() < sample_size / max_sequences:
            sample_sequences.append(seq)
            
        for aa in seq:
            if aa in STANDARD_AMINO_ACIDS:
                aa_counts[aa] += 1
    
    # Force garbage collection
    gc.collect()
    
    return aa_counts, total_length, sample_sequences, sequence_lengths

def bootstrap_iteration(args):
    """Perform a single bootstrap iteration with reduced memory usage."""
    sequences, func, aa, _ = args
    
    # Use reservoir sampling for large datasets
    sample_size = min(1000, len(sequences))
    indices = np.random.choice(len(sequences), size=sample_size, replace=True)
    resampled = [sequences[i] for i in indices]
    
    if aa is not None:
        return func(resampled, aa)
    else:
        return func(resampled)

def bootstrap_error(sequences, func, aa=None, n_bootstrap=2):
    """Estimate error using bootstrapping with reduced iterations."""
    if not sequences:
        return 0
        
    # Limit sequence sample size for bootstrapping
    max_sample = min(5000, len(sequences))
    if len(sequences) > max_sample:
        sequences = random.sample(sequences, max_sample)
    
    # Prepare arguments for bootstrap_iteration
    args = [(sequences, func, aa, i) for i in range(n_bootstrap)]
    
    # Use fewer processes for bootstrapping
    n_processes = min(n_bootstrap, max(1, cpu_count() // 2))
    with Pool(n_processes) as pool:
        values = pool.map(bootstrap_iteration, args)
    
    return np.std(values) if values else 0

def process_group_parallel(directory, group_name):
    """Process a group of FASTA files in parallel with resource constraints."""
    aa_counts = defaultdict(int)
    total_length = 0
    all_sequences = []
    all_lengths = []
    
    print(f"Processing directory: {directory} for group: {group_name}")
    file_paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if (file.startswith(group_name) or group_name == "Selected Organisms") and (file.endswith('.fasta.gz') or file.endswith('.fasta')):
                file_paths.append(os.path.join(root, file))
    
    # Limit the number of parallel processes
    max_sequences = 5000  # Limit sequences per file
    sample_size = 500  # Sample size for bootstrapping
    n_processes = max(1, min(cpu_count() // 2, 4))  # Use at most 4 processes or half CPU cores
    
    print(f"Using {n_processes} processes for parallel processing")
    
    # Prepare arguments for process_file
    args = [(file_path, max_sequences, sample_size) for file_path in file_paths]
    
    # Use multiprocessing with limited processes
    with Pool(n_processes) as pool:
        results = pool.map(process_file, args)
    
    for counts, length, sequences, lengths in results:
        for aa, count in counts.items():
            aa_counts[aa] += count
        total_length += length
        all_sequences.extend(sequences)
        all_lengths.extend(lengths)
        
        # Clear memory after each batch
        gc.collect()
    
    if not all_sequences and not all_lengths:
        print(f"No sequences found for {group_name} in {directory}")
        return {}, {}, 0, 0
    
    # Calculate average aa content
    aa_content = {}
    if total_length > 0:
        aa_content = {aa: (count / total_length) * 100 for aa, count in aa_counts.items()}
    
    print(f"Calculating amino acid content errors for {group_name}...")
    aa_content_error = {}
    for aa in STANDARD_AMINO_ACIDS:
        aa_content_error[aa] = bootstrap_error(all_sequences, calculate_aa_percentage, aa)
        # Clear memory after each amino acid
        gc.collect()
    
    # Calculate average protein length
    avg_length = np.mean(all_lengths) if all_lengths else 0
    print(f"Calculating protein length error for {group_name}...")
    length_error = bootstrap_error(all_lengths, np.mean, None)
    
    # Clear memory before returning
    all_sequences = []
    all_lengths = []
    gc.collect()
    
    return aa_content, aa_content_error, avg_length, length_error

def main():
    """Main function to process data and generate tables."""
    # Set lower memory usage target
    import sys
    if sys.platform.startswith('linux'):
        try:
            import resource
            # Set soft limit to 4GB
            resource.setrlimit(resource.RLIMIT_AS, (4 * 1024 * 1024 * 1024, -1))
        except (ImportError, resource.error):
            pass
    
    base_dir = os.path.abspath(os.path.dirname(__file__))
    selected_organisms_dir = os.path.join(base_dir, "selected_organisms")
    random_selected_dir = os.path.join(base_dir, "random_selected")
    
    print(f"Base directory: {base_dir}")
    print(f"Selected organisms directory: {selected_organisms_dir}")
    print(f"Random selected directory: {random_selected_dir}")
    
    # Process selected organisms (a)
    print("\nProcessing selected organisms (a)...")
    aa_content_a, aa_content_error_a, avg_length_a, length_error_a = process_group_parallel(selected_organisms_dir, "Selected Organisms")
    
    # Force garbage collection between major operations
    gc.collect()
    
    # Process randomly selected kingdoms (c)
    print("\nProcessing randomly selected kingdoms (c)...")
    kingdoms = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
    aa_content_c = {}
    aa_content_error_c = {}
    avg_length_c = {}
    length_error_c = {}
    
    for kingdom in kingdoms:
        print(f"\nProcessing {kingdom} files in directory: {random_selected_dir}")
        aa_content_c[kingdom], aa_content_error_c[kingdom], avg_length_c[kingdom], length_error_c[kingdom] = process_group_parallel(random_selected_dir, kingdom)
        # Force garbage collection after each kingdom
        gc.collect()
    
    # Create tables using PrettyTable
    table1 = PrettyTable()
    table1.field_names = ["Amino Acid", "Average Content (%)", "Error"]
    for aa in sorted(STANDARD_AMINO_ACIDS):
        table1.add_row([aa, f"{aa_content_a.get(aa, 0):.2f}", f"{aa_content_error_a.get(aa, 0):.2f}"])
    
    table2 = PrettyTable()
    table2.field_names = ["Group", "Average Length", "Error"]
    table2.add_row(["Selected Organisms", f"{avg_length_a:.2f}", f"{length_error_a:.2f}"])
    
    table3 = PrettyTable()
    table3.field_names = ["Kingdom", "Amino Acid", "Average Content (%)", "Error"]
    for kingdom in kingdoms:
        for aa in sorted(STANDARD_AMINO_ACIDS):
            table3.add_row([kingdom, aa, f"{aa_content_c[kingdom].get(aa, 0):.2f}", f"{aa_content_error_c[kingdom].get(aa, 0):.2f}"])
    
    table4 = PrettyTable()
    table4.field_names = ["Kingdom", "Average Length", "Error"]
    for kingdom in kingdoms:
        table4.add_row([kingdom, f"{avg_length_c[kingdom]:.2f}", f"{length_error_c[kingdom]:.2f}"])
    
    print("\nTable 1: Average Amino Acid Content for Selected Organisms (a)")
    print(table1)
    
    print("\nTable 2: Average Protein Length for Selected Organisms (a)")
    print(table2)
    
    print("\nTable 3: Average Amino Acid Content for Randomly Selected Kingdoms (c)")
    print(table3)
    
    print("\nTable 4: Average Protein Length for Randomly Selected Kingdoms (c)")
    print(table4)

if __name__ == '__main__':
    main()