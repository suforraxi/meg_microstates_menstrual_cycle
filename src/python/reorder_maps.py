import glob
import os
import numpy as np
import pandas as pd
from micros_based_stats.examples.plot_brain_aal import plot_brains
from statsmodels.stats.multitest import multipletests # Import for FDR correction
import seaborn as sns

import matplotlib.pyplot as plt

# Apply FDR correction per N_Microstate
def apply_fdr_correction(group, alpha=0.05):
    """
    Apply False Discovery Rate (FDR) correction to a group of p-values.

    This function adjusts the p-values in the input group using the Benjamini-Hochberg
    procedure to control the false discovery rate. It adds two new columns to the input
    DataFrame: 'FDR_p', which contains the FDR-corrected p-values, and 'Significant_FDR',
    which indicates whether each p-value is significant after the correction.

    Parameters:
    -----------
    group : pandas.DataFrame
        A DataFrame containing a column named "p-value" with the p-values to be corrected.
    alpha : float, optional
        The significance level to control the false discovery rate (default is 0.05).

    Returns:
    --------
    pandas.DataFrame
        The input DataFrame with two additional columns:
        - 'FDR_p': The FDR-corrected p-values.
        - 'Significant_FDR': A boolean column indicating whether each p-value is
            significant after FDR correction.

    Notes:
    ------
    This function uses the `multipletests` function from the `statsmodels.stats.multitest`
    module with the 'fdr_bh' method for FDR correction.

    Example:
    --------
    >>> import pandas as pd
    >>> from statsmodels.stats.multitest import multipletests
    >>> data = pd.DataFrame({"p-value": [0.01, 0.04, 0.03, 0.2, 0.5]})
    >>> corrected_data = apply_fdr_correction(data)
    >>> print(corrected_data)
    """
    _, fdr_corrected_pvals, _, _ = multipletests(
        group["p-value"], alpha=alpha, method="fdr_bh"
    )
    group["FDR_p"] = fdr_corrected_pvals
    group["Significant_FDR"] = group["FDR_p"] < alpha
    return group


def reorder_maps(in_folder, base_output_dir):
    all_files = sorted(
        glob.glob(
            os.path.join(
                in_folder,
                "**_microstates_segmentation",
                "microstates_all_subjects_results.npz",
            )
        )
    )
    max_F_x_k = {}
    for f in all_files:
        n_microstates = os.path.basename(os.path.dirname(f)).split("_")[0]
        # 1. Use a wildcard '*' to get all folders matching the pattern
        maps_file = os.path.join(
            in_folder,
            f"{n_microstates}_microstates_segmentation",
            "microstates_all_subjects_results.npz",
        )

        # Load the .npz file
        data = np.load(maps_file, allow_pickle=True)
        # Extract the microstate_maps matrix
        if "Microstates" in data:
            microstate_maps = data["Microstates"]
            print(f"Loaded microstate_maps with shape: {microstate_maps.shape}")
        else:
            raise KeyError("The key 'Microstates' was not found in the .npz file.")

        # load the anova results
        anova_table_f = os.path.join(
            in_folder,
            f"{n_microstates}_backfitted_microstates",
            "combined",
            "anova_results.csv",
        )
        df_anova = pd.read_csv(anova_table_f)

        # Apply FDR correction to the p-values
        df_anova = apply_fdr_correction(df_anova, alpha=0.05)

        df_anova = df_anova.sort_values(by="F", ascending=False)

        print(df_anova["State"])

        # Extract the order of states from df_anova
        state_order = df_anova['State'].values
        max_F_x_k[int(n_microstates)] = df_anova['F'].max()
        # Reorder the rows of microstate_maps based on state_order
        microstate_maps_ordered = microstate_maps[state_order, :]

        # Save the reordered microstate_maps to a new file
        output_file = maps_file.replace(
            "microstates_all_subjects_results.npz",
            "microstates_all_subjects_results_ordered.npz",
        )

        np.savez(output_file, Microstates=microstate_maps_ordered)

        print(f"Reordered microstate_maps saved to {output_file}")

        # Plot the reordered microstate maps
        #plot_brains(
        #    output_file,
        #    os.path.join(base_output_dir, f"{n_microstates}_microstates_ordered"),
        #)
    df = pd.DataFrame(list(max_F_x_k.items()), columns=['K', 'Max_F'])
    sorted_df = df.sort_values(by='K')
    sorted_df.to_csv(os.path.join(base_output_dir, "max_F_x_k.csv"), index=False)
def corr_inv_dot(source_vector, target_map):
    """
    Compute polarity-invariant correlation using dot product between a source vector and rows of target_map.

    Parameters:
    -----------
    source_vector : np.ndarray
        A 1D array representing the source vector (shape: (x,)).
    target_map : np.ndarray
        A 2D array where each row is a target vector (shape: (y, x)).

    Returns:
    --------
    np.ndarray
        A 1D array of shape (y,) containing the maximum correlation values for the source vector
        with each target row, considering polarity invariance.
    """
    # Normalize source_vector and rows of target_map
    source_norm = source_vector / np.linalg.norm(source_vector)
    target_norm = target_map / np.linalg.norm(target_map, axis=1, keepdims=True)

    # Compute dot product for original and inverted polarity
    corr = np.dot(target_norm, source_norm)
    corr_inv = np.dot(target_norm, -source_norm)

    # Take the maximum correlation for polarity invariance
    return np.maximum(corr, corr_inv)


def correlate_maps(maps_file_source,
                   source_rois,
                   in_folder,
                   base_output_dir):
     # Load the .npz file
    data = np.load(maps_file_source, allow_pickle=True)
    # Extract the microstate_maps matrix
    if "Microstates" in data:
        source_maps = data["Microstates"]
        print(f"Loaded microstate_maps with shape: {source_maps.shape}")
    else:
        raise KeyError("The key 'Microstates' was not found in the .npz file.")

    source_matrix = source_maps[source_rois, :]

    all_files = sorted(
        glob.glob(
            os.path.join(
                in_folder,
                "**_microstates_segmentation",
                "microstates_all_subjects_results_ordered.npz",
            )
        )
    )
    all_corr_m = []
    sig_k = []
    n_sig_x_k = {}
    for f in all_files:
        n_microstates = os.path.basename(os.path.dirname(f)).split("_")[0]
        data = np.load(f, allow_pickle=True)
        target_map = None
        # Extract the microstate_maps matrix
        if "Microstates" in data:
            target_map = data["Microstates"]
            print(f"Loaded microstate_maps with shape: {target_map.shape}")
        else:
            raise KeyError("The key 'Microstates' was not found in the .npz file.")

        anova_table_f = os.path.join(
                in_folder,
                f"{n_microstates}_backfitted_microstates",
                "combined",
                "anova_results.csv",
            )
        df_anova = pd.read_csv(anova_table_f)

        # Apply FDR correction to the p-values
        df_anova = apply_fdr_correction(df_anova, alpha=0.05)

        df_anova = df_anova.sort_values(by="F", ascending=False)
        # Rename states based on their order
        df_anova['State'] = range(len(df_anova))
        
        df_anova.to_csv(anova_table_f.replace(
            "anova_results.csv",
            "anova_results_ordered.csv"
            ), index=False
        )

        n_sig_states = df_anova['Significant_FDR'].sum()
        if n_sig_states > 0:
            n_sig_x_k[int(n_microstates)] = n_sig_states
            sig_k.append(int(n_microstates))
            print(f"Number of significant states after FDR correction: {n_sig_states}")
            #Compute the correlation between source_maps and target_map
            corr_m = np.zeros((source_matrix.shape[0], target_map.shape[0]))
            for i in range(source_matrix.shape[0]):
                source_vector = source_matrix[i, :]
                corr_m[i, :] = corr_inv_dot(source_vector, target_map)
            
            all_corr_m.append(corr_m)

            fig = plt.figure(figsize=(10, 8))
            plt.imshow(corr_m, aspect='auto', cmap='viridis')
            plt.colorbar(label='Spatial Correlation')
            plt.xticks(ticks=np.arange(target_map.shape[0]), 
                    labels=[f'HDM {i}' for i in range(target_map.shape[0])],
                    rotation=90,
                    fontsize=10
            )
            plt.yticks(ticks=np.arange(source_matrix.shape[0]), 
                    labels=[f'HDM {i}' for i in range(source_matrix.shape[0])],
                    fontsize=10
            )
            plt.xlabel('Target Maps', fontsize=12)
            plt.ylabel('Source Maps (k=13 optimal solution)', fontsize=12)
            plt.title(f'K = {n_microstates} clusters\nNumber of significant microstates across MC phases (FDR corrected): {n_sig_states}'
            )

            # Create a red empty box for the boundary of significant columns
            for col in range(n_sig_states):
                plt.gca().add_patch(plt.Rectangle(
                    (col - 0.5, -0.5), 1, source_matrix.shape[0],
                    edgecolor='red', fill=False, lw=4, joinstyle='miter'  # Uniform edge thickness
                ))          

            plt.savefig(os.path.join(base_output_dir, 
                                    f'{n_microstates}_microstates_correlation.png'
                                    )
            )
            plt.close(fig)
    
    max_n_hdm = max([m.shape[1] for m in all_corr_m])
    m = np.zeros((max_n_hdm, len(all_corr_m)))
    
    source_map = [0, 1] # Assuming you want to use the first two source maps for correlation
    for source_idx in source_map:
        for i in range(len(all_corr_m)):
                n_hdm = all_corr_m[i].shape[1]
                m[:n_hdm, i] = all_corr_m[i][source_idx, :n_hdm]
        
        heat_m_df = pd.DataFrame(m, columns=[f'K={k}' for k in sig_k])
        # Reorder the columns in ascending order according to k
        heat_m_df = heat_m_df[sorted(heat_m_df.columns, key=lambda col: int(col.split('=')[1]))]

        # Set the index names as 'HDM 0', 'HDM 1', etc.
        heat_m_df.index = [f'HDM {i}' for i in range(heat_m_df.shape[0])]

        # Create a mask for zeros to make them black and hide the numbers
        mask = (heat_m_df == 0)

        plt.figure(figsize=(12, 10))  # Adjust the figure size as needed
        hm = sns.heatmap(heat_m_df,
                    cmap='viridis',
                    annot=True,
                    annot_kws={"size": 12,
                               "weight": 'bold'
                        },
                    fmt=".2f",
                    mask=mask,
                    cbar_kws={'label': 'Spatial Correlation',
                              'ticks': []
                              }
                )

        # Overlay black boxes for zeros
        for (i, j), val in np.ndenumerate(heat_m_df.values):
            if val == 0:
                plt.gca().add_patch(plt.Rectangle((j, i), 1, 1, color='white'))  # Add white box to hide the number
                plt.gca().add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=0.5))  # Add visible borders

        # Add borders to the bottom of the heatmap
        plt.gca().add_patch(plt.Rectangle((-0.5, heat_m_df.shape[0] - 0.1), heat_m_df.shape[1], 1, 
                                          fill=False, edgecolor='black', lw=1))

        plt.title("")
        plt.xlabel("K solutions", fontsize=20)
        plt.ylabel("HDM maps", fontsize=20)
        plt.xticks(rotation=45, fontsize=18)
        plt.yticks(rotation=0, fontsize=18)
        hm.figure.axes[-1].yaxis.label.set_size(18) 
        #plt.show()
        save_path = os.path.join(base_output_dir, f"{source_idx}_heat_map_matrices")
        plt.savefig(save_path)
        plt.close()
        
        # Create a new DataFrame for n_sig_x_k values
        n_sig_x_k_df = pd.DataFrame.from_dict(n_sig_x_k, orient='index',
                                              columns=['n_sig_states'])
        n_sig_x_k_df.index.name = 'K'

        # Order the rows based on K
        n_sig_x_k_df = n_sig_x_k_df.sort_index()
        #n_sig_x_k_df['K-clusters'] = n_sig_x_k_df.index  # Convert K to string for better display
        #col_order = ['K-clusters','n_sig_states']
        #n_sig_x_k_df = n_sig_x_k_df[col_order]
        
        # Update the heatmap for n_sig_x_k values
        plt.figure(figsize=(6, 16))  # Keep the compact width
        sns.heatmap(n_sig_x_k_df,
                    cmap=sns.color_palette("coolwarm", as_cmap=True),
                    annot=True,
                    annot_kws={"size": 25, "weight": 'bold'},  # Adjust annotation font size and weight
                    fmt="d",
                    linewidths=0.5,
                    linecolor='black',
                    cbar=False,
                    cbar_kws={'label': 'Number of Significant HDMs\nafter FDR Correction', 
                              'ticks': []}
                              )  # Remove ticks
        plt.title("",
                  fontsize=18,
                  fontweight='bold',
                  pad=10
                  )
        plt.xlabel("")  # Remove column labels
        plt.xticks([])  # Hide x-axis ticks
        plt.yticks(rotation=0, 
                   fontsize=20,
                   fontweight='bold')  # Rotate y-ticks to horizontal and adjust font size
        plt.ylabel("",
                   rotation=0,
                   fontsize=22,
                   fontweight='bold',
                   loc='top'
                   )

        # Save the heatmap
        save_path = os.path.join(base_output_dir, "n_sig_x_k_heatmap.png")
        plt.savefig(save_path)
        plt.close()
        n_sig_x_k_df.to_csv(os.path.join(base_output_dir,
                                         "n_sig_x_k.csv"),
                            index=True
                            )



# Path to the .npz file
in_folder = os.path.join(os.getcwd(), "results", "24_1000")
base_output_dir = os.path.join(os.getcwd(), "results", "24_1000", "ordered_brains")
corr_output_dir = os.path.join(os.getcwd(), "results", "24_1000", "correlation_matrices")

if not os.path.exists(base_output_dir):
    os.makedirs(base_output_dir)
if not os.path.exists(corr_output_dir):
    os.makedirs(corr_output_dir)

max_F_x_k = reorder_maps(in_folder, base_output_dir)
n_microstates = 13
maps_file_source = os.path.join(
    in_folder,
    f"{n_microstates}_microstates_segmentation",
    "microstates_all_subjects_results_ordered.npz",
)
correlate_maps(maps_file_source,
                source_rois=[0, 1],
                in_folder=in_folder,
                base_output_dir=corr_output_dir
                )

