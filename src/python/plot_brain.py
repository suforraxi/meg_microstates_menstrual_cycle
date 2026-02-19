import glob 

import os
from micros_based_stats.examples.plot_brain_aal import plot_brains

# Path to the .npz file
in_folder = os.path.join(os.getcwd(), 'results', '24_1000')
base_output_dir = os.path.join(os.getcwd(), 'src', 'r', 'figures', 'brains')

if not os.path.exists(base_output_dir):
    os.makedirs(base_output_dir)

# 1. Use a wildcard '*' to get all folders matching the pattern
all_files = glob.glob(os.path.join(in_folder, '*_microstates_segmentation', 'microstates_all_subjects_results.npz'))

# 2. Filter for the numeric range 13
filtered_files = [
    f for f in all_files 
    if 13 == int(os.path.basename(os.path.dirname(f)).split('_')[0]) 
]
for f in filtered_files:
    n_microstates = os.path.basename(os.path.dirname(f)).split('_')[0]
    out_dir = os.path.join(base_output_dir, f'{n_microstates}_microstates')
    plot_brains(f, out_dir)
