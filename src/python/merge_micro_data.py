import os
import pandas as pd

r_folder = './src/r/data/'

# Load the two tables
file_0 = os.path.join(r_folder, "data_microstate_7_for_R.csv")
file_1 = os.path.join(r_folder, "data_microstate_2_for_R.csv")

df_0 = pd.read_csv(file_0)
df_1 = pd.read_csv(file_1)

# Rename 'Microstate_z' and 'Occurrences' to distinguish between the two files
df_0 = df_0.rename(columns={'Microstate_z': 'Microstate_z_7', 'Occurrences': 'Occurrences_7'})
df_1 = df_1.rename(columns={'Microstate_z': 'Microstate_z_2', 'Occurrences': 'Occurrences_2'})

# Merge the tables on all common columns, keeping all columns and adding the renamed ones
merged_df = pd.merge(
    df_0,
    df_1,
    on=[col for col in df_0.columns if col in df_1.columns and col not in ['Microstate_z_7', 'Occurrences_7', 'Microstate_z_2', 'Occurrences_2']],
    how='outer',  # Use outer join to keep all rows and columns
    suffixes=('_7', '_2')  # Ensure unique suffixes for any other overlapping columns
)

# Save the merged table to a new CSV file
merged_output_file = os.path.join(r_folder, "merged_data_microstates_for_R.csv")
merged_df.to_csv(merged_output_file, index=False)
print(f"Merged DataFrame saved to: {merged_output_file}")