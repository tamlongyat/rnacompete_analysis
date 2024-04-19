import os

import matplotlib.pyplot as plt
import numpy as np

def plot_scatter(kmer_zscore_df, scatter_folder):
    """
    Plot and save the scatter plot for the k-mer Z-score correlation between set A and B
    for a list of proteins.
    
    Parameters
    ----------
    kmer_zscore_df : pd.DataFrame
        DataFrame containing the k-mer Z-scores for different proteins and sets.
    scatter_folder : str
        Output folder containing the scatter plots.
    
    """
    
    # Iterate over the proteins
    for hybid in [e.split('_')[0] for e in kmer_zscore_df.columns[::3]]:
    
        # Get the Z-scores
        zscore_a_list = kmer_zscore_df[f'{hybid}_a'].to_list()
        zscore_b_list = kmer_zscore_df[f'{hybid}_b'].to_list()
        
        # Generate the plot
        plt.figure(figsize=(4, 4))
        plt.scatter(zscore_a_list, zscore_b_list, color='black', s=5, zorder=3)
        
        # Set the limits
        lim_min = min(zscore_a_list + zscore_b_list) - 1
        lim_max = max(zscore_a_list + zscore_b_list)
        lim_max = np.ceil(lim_max / 5) * 5
        plt.xlim(lim_min, lim_max)
        plt.ylim(lim_min, lim_max)
        
        # Set the ticks
        plt.xticks(np.arange(0, lim_max + 5, 5))
        plt.yticks(np.arange(0, lim_max + 5, 5))
        
        # Set the labels
        plt.xlabel('Set A', weight='bold')
        plt.ylabel('Set B', weight='bold')
        plt.grid(color='lightgray', zorder=0)
        
        # Save the figure
        scatter_path = os.path.join(scatter_folder, f'{hybid}.png')
        plt.savefig(scatter_path)
        plt.close()
