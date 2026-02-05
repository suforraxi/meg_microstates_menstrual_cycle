from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import zscore


def process_microstate_psi(microstate_file, psi_file, output_folder):
    """
    Processes microstate occurrences and psychological test data, creating a combined DataFrame for each Microstate.

    Parameters:
        microstate_file (str): Path to the microstate occurrences file.
        psi_file (str): Path to the psychological test data file.
        output_folder (str): Folder to save the output files.
    """
    # Load the data
    microstate_df = pd.read_csv(microstate_file)
    psi_df = pd.read_csv(psi_file)

    # Reshape the psychological test data to long format
    psi_long = pd.melt(
        psi_df,
        id_vars=['cod', 'sub'],
        value_vars=[
            'BAI_S1', 'ROSENBERG_S1', 'BDI_S1', 'WELLBEING_S1', 'Autonomy_S1', 'EnvironmentalMastery_S1',
            'PersonalGrowth_S1', 'PositiveRelationswithOthers_S1', 'PurposeinLife_S1', 'Self-Acceptance_S1',
            'BAI_S2', 'ROSENBERG_S2', 'BDI_S2', 'WELLBEING_S2', 'Autonomy_S2', 'EnvironmentalMastery_S2',
            'PersonalGrowth_S2', 'PositiveRelationswithOthers_S2', 'PurposeinLife_S2', 'Self-Acceptance_S2',
            'BAI_S3', 'ROSENBERG_S3', 'BDI_S3', 'WELLBEING_S3', 'Autonomy_S3', 'EnvironmentalMastery_S3',
            'PersonalGrowth_S3', 'PositiveRelationswithOthers_S3', 'PurposeinLife_S3', 'Self-Acceptance_S3'
        ],
        var_name='Psi_Session',
        value_name='Psi_Value'
    )

    # Extract the psychological test label and session number
    psi_long['Psi_Label'] = psi_long['Psi_Session'].str.split('_').str[0]
    psi_long['Session'] = psi_long['Psi_Session'].str.split('_').str[1]

    psi_long['Session'] = psi_long['Session'].map({'S1': 1, 'S2': 2, 'S3': 3})
    # Drop the Psi_Session column as it's no longer needed
    psi_long = psi_long.drop(columns=['Psi_Session'])

    # Iterate over each Microstate and create a separate DataFrame
    unique_microstates = microstate_df['Microstate'].unique()

    for microstate in unique_microstates:
        # Filter the microstate data for the current Microstate
        microstate_filtered = microstate_df[microstate_df['Microstate'] == microstate]
        
        # Merge the filtered microstate data with the reshaped psychological test data
        merged_df = pd.merge(
            microstate_filtered,
            psi_long,
            on=['sub', 'Session'],
            how='inner'
        )
        
        # Select and reorder the desired columns
        final_df = merged_df[['sub', 'Session', 'Psi_Label', 'Psi_Value', 'Occurrences_zscore', 'Occurrences']]
        
        # Save the final DataFrame to a CSV file
        output_file = os.path.join(output_folder, f'combined_psi_microstates_Microstate_{microstate}.csv')
        final_df.to_csv(output_file, index=False)
        print(f"Final DataFrame for Microstate {microstate} created and saved to:", output_file)

def process_microstate_hormones(n_peaks, n_microstates, microstate_file, hormones_file, output_folder):
    """
    Processes microstate occurrences and hormone data, creating a combined DataFrame for each Microstate.

    Parameters:
        n_peaks (int): Number of peaks per acquisition.
        n_microstates (int): Number of microstates.
        microstate_file (str): Path to the microstate occurrences file.
        hormones_file (str): Path to the hormones data file.
        output_folder (str): Folder to save the output files.
    """
    # Load the data
    microstate_df = pd.read_csv(microstate_file)
    hormones_df = pd.read_csv(hormones_file, sep=',')

    # Reshape the hormones data to long format
    hormones_long = pd.melt(
        hormones_df,
        id_vars=['cod', 'sub'],
        value_vars=['P_S1', 'P_S2', 'P_S3', 'LH_S1', 'LH_S2', 'LH_S3', 'FSH_S1', 'FSH_S2', 'FSH_S3', 'E_S1', 'E_S2', 'E_S3'],
        var_name='Hormone_Session',
        value_name='Hormone_Value'
    )

    # Extract Hormone and Session from the Hormone_Session column
    hormones_long['Hormone'] = hormones_long['Hormone_Session'].str.split('_').str[0]
    hormones_long['Session'] = hormones_long['Hormone_Session'].str.split('_').str[1]

    # Map Session values from 'S1', 'S2', 'S3' to '1', '2', '3'
    hormones_long['Session'] = hormones_long['Session'].map({'S1': 1, 'S2': 2, 'S3': 3})

    # Drop the Hormone_Session column as it's no longer needed
    hormones_long = hormones_long.drop(columns=['Hormone_Session'])

    # Iterate over each Microstate and create a separate DataFrame
    unique_microstates = microstate_df['Microstate'].unique()

    for microstate in unique_microstates:
        # Filter the microstate data for the current Microstate
        microstate_filtered = microstate_df[microstate_df['Microstate'] == microstate]
        
        # Merge the filtered microstate data with the reshaped hormone data
        merged_df = pd.merge(
            microstate_filtered,
            hormones_long,
            on=['sub', 'Session'],
            how='inner'
        )
        
        # Select and reorder the desired columns
        final_df = merged_df[['sub', 'Session', 'Hormone', 'Hormone_Value', 'Occurrences','Occurrences_zscore']]
        
        # Save the final DataFrame to a CSV file
        output_file = os.path.join(output_folder, f'combined_hormones_microstates_Microstate_{microstate}.csv')
        final_df.to_csv(output_file, index=False)
        print(f"Final DataFrame for Microstate {microstate} created and saved to:", output_file)

# Function to perform ML-PCA
def ml_pca(data, group_col, components=1):
    """
    Perform Multilevel PCA (ML-PCA) to separate between-subject and within-subject variance.

    Parameters:
        data (pd.DataFrame): DataFrame with repeated measures data.
        group_col (str): Column name for the grouping variable (e.g., 'sub').
        components (int): Number of principal components to extract.

    Returns:
        dict: Dictionary with within-subject PCA results.
    """
    # Center the data within each group (within-subject variance)
    data_centered = data.groupby(group_col).transform(lambda x: x - x.mean())

    # Z-score the data (standardize to mean=0, std=1)
    data_zscored = data_centered.apply(zscore, axis=0)

    # Compute within-subject PCA
    scaler_within = StandardScaler()
    data_within_scaled = scaler_within.fit_transform(data_zscored)
    pca_within = PCA(n_components=components)
    within_pcs = pca_within.fit_transform(data_within_scaled)

    return {
        "within": {"pca": pca_within, "pcs": within_pcs, "explained_variance": pca_within.explained_variance_ratio_},
    }

def plot_occurrences_per_subject(df, output_folder):
    """
    Plots Occurrences for each subject across different sessions, with lines connecting the values,
    and saves the plot to the specified folder. Adds a black dotted horizontal line at zero and ensures
    x-ticks are properly aligned with the plotted points.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the relevant columns.
        output_folder (str): The folder where the plot will be saved.
    """
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Ensure Session is treated as a categorical variable
    df['Session'] = df['Session'].astype('category')

    # Sort the DataFrame by subject to ensure proper alignment
    df = df.sort_values(by='sub')

    # Create the plot
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=df,
        x='sub',  # Subject on the x-axis
        y='Occurrences',  # Occurrences on the y-axis
        hue='Session',  # Different lines for each session
        marker='o',  # Markers for each point
        palette='tab10'
    )

    # Add a black dotted horizontal line at zero
    plt.axhline(0, color='black', linestyle='dotted', linewidth=1)

    # Customize the plot
    plt.title('Occurrences per Subject Across Sessions', fontsize=16)
    plt.xlabel('Subject', fontsize=14)
    plt.ylabel('Occurrences', fontsize=14)

    # Set x-ticks to match the unique subjects
    #unique_subjects = df['sub'].unique()
    #plt.xticks(ticks=range(len(unique_subjects)), labels=unique_subjects, rotation=45, fontsize=10)

    plt.legend(title='Session', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.tight_layout()

    # Save the plot
    output_file = os.path.join(output_folder, "occurrences_per_subject.png")
    plt.savefig(output_file, bbox_inches='tight')
    print(f"Occurrences plot saved to: {output_file}")

    # Show the plot (optional, can be commented out if not needed)
    plt.close()

# Example usage
n_peaks = '1000'
n_microstates = 13
base_folder= os.path.join(os.getcwd(), 'results', n_peaks)
data_folder = os.path.join(base_folder, f'{n_microstates}_backfitted_microstates', 'combined')
microstates = [2, 7]
psi_labels = ['Autonomy', 'EnvironmentalMastery',
              'PersonalGrowth', 'PositiveRelationswithOthers', 
              'PurposeinLife', 'Self-Acceptance']
hormones = ['E', 'P', 'FSH','LH']

microstate_file = os.path.join(data_folder, "aggregated_microstate_occurrences.csv")
hormones_file = os.path.join(os.getcwd(), 'data', "hormones_data.csv")
psi_file = os.path.join(os.getcwd(), 'data', "psi_data.csv")
output_folder = data_folder
r_folder = os.path.join(os.getcwd(), "src", "r", "data")
process_microstate_hormones(n_peaks, n_microstates, microstate_file, hormones_file, output_folder)
process_microstate_psi(microstate_file, psi_file, output_folder)

microstate_df = pd.read_csv(microstate_file)
microstate_df = microstate_df[microstate_df['Microstate'].isin(microstates)]
plot_occurrences_per_subject(microstate_df, output_folder)

for microstate in microstates:
    # Load the data
    hormones_file = os.path.join(data_folder, f'combined_hormones_microstates_Microstate_{microstate}.csv')
    psi_file = os.path.join(data_folder, f'combined_psi_microstates_Microstate_{microstate}.csv')
    
    hormones_df = pd.read_csv(hormones_file)
    psi_df = pd.read_csv(psi_file)

    # Merge the hormones and psychological data on 'sub' and 'Session'
    merged_df = pd.merge(hormones_df, psi_df, on=['sub', 'Session'], how='inner', suffixes=('', '_psi'))
    merged_df = merged_df.loc[:, ~merged_df.columns.str.endswith('_psi')]

    # Pivot psychological data
    pivot_psych = merged_df.pivot_table(
        index=['sub', 'Session', 'Occurrences_zscore', 'Occurrences'],
        columns='Psi_Label',
        values='Psi_Value'
    ).reset_index()

    # Pivot the hormone data to create separate columns for each hormone
    pivot_hormones = merged_df.pivot_table(
        index=['sub', 'Session'],
        columns='Hormone',
        values='Hormone_Value'
    ).reset_index()

    filtered_data = pd.merge(pivot_hormones, pivot_psych, on=['sub', 'Session'], how='left')

    df = filtered_data.drop_duplicates().copy()

    # Perform ML-PCA on hormones
    hormone_data = df[['sub', 'Session'] + hormones]
    
    ml_pca_hormones = ml_pca(hormone_data, group_col='sub', components=1)

    # Add within-subject PC1 for hormones directly
    df['Hormones_PC1_within'] = ml_pca_hormones["within"]["pcs"][:, 0]

    # Print variance explained by within PC1 for hormones
    print(f"Variance explained by within-subject PC1 for hormones: {ml_pca_hormones['within']['explained_variance'][0]:.2%}")

    # Perform ML-PCA on psychological scores
    psych_data = df[['sub', 'Session'] + psi_labels]
    
    ml_pca_psy = ml_pca(psych_data , group_col='sub', components=1)

    # Add within-subject PC1 for hormones directly
    df['Psy_PC1_within'] =  ml_pca_psy["within"]["pcs"][:, 0]

    # Print variance explained by within PC1 for hormones
    print(f"Variance explained by within-subject PC1 for psy: { ml_pca_psy['within']['explained_variance'][0]:.2%}")


    # Ensure Session is treated as a categorical variable
    df['Session'] = df['Session'].astype('category')

    #df['WELLBEING'] = psych_data['WELLBEING']

    # Filter the psychological data to include only 'WELLBEING'
    #psych_data = df[['sub', 'Session'] + psi_labels]

    # Use raw values directly 'WELLBEING'
    #df['WELLBEING'] = psych_data['WELLBEING']

    # Save the DataFrame to a CSV file for debugging or further analysis
    df_output_file = os.path.join(output_folder, f"data_microstate_{microstate}.csv")
    df.to_csv(df_output_file, index=False)
    print(f"DataFrame saved to: {df_output_file}")
    df_output_file_r = os.path.join(r_folder, f"data_microstate_{microstate}_for_R.csv")
    df = df.rename(columns={'Occurrences_zscore': 'Microstate_z'})
    df.to_csv(df_output_file_r, index=False)






