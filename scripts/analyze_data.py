import os
import sys

from absl import logging
import numpy as np
import pandas as pd

sys.path.append('subscripts')
from compute_probe_zscore import compute_probe_zscore
from compute_kmer_zscore import compute_kmer_zscore
from plot_scatter import plot_scatter
from align_top_kmer import align_top_kmer
from compute_pwm import compute_pwm
from plot_logo import plot_logo
from run_classifier import run_classifier


def get_top_kmer(kmer_zscore_df):
    """
    Get the top k-mer sand Z-scores for all proteins and sets.
    
    Parameters
    ----------
    kmer_zscore_df : pd.DataFrame
        DataFrame containing the k-mer Z-scores for different proteins and sets.
    
    Returns
    -------
    top_kmer_mat : np.ndarray
        Matrix of the top 10 k-mers for the proteins in different sets.
    top_zscore_mat : np.ndarray
        Matrix of the Z-scores of the top 10 k-mers for the proteins in different sets.
    """
    
    # Iterate over the proteins and sets
    top_kmer_mat = []
    top_zscore_mat = []
    for hybid_set in kmer_zscore_df.columns:
        
        # Get the top k-mers and Z-scores
        tmp_df = kmer_zscore_df[hybid_set].sort_values(ascending=False)[:10]
        top_kmer_mat.append(tmp_df.index.to_list())
        top_zscore_mat.append(tmp_df.to_list())
    top_kmer_mat = np.array(top_kmer_mat)
    top_zscore_mat = np.array(top_zscore_mat)
    return top_kmer_mat, top_zscore_mat


def get_summary(batch_metadata_df, top_kmer_mat, top_zscore_mat, aligned_top_kmer_mat, class_list, prob_list):
    
    # Get the column names
    col_list = ['rncmpt', 'gene_name', 'species', 'class', 'prob']
    for prefix in ['kmer', 'zscore', 'aligned_kmer']:
        for set_name in ['a', 'b', 'ab']:
            for idx in range(1, 11):
                col_list += [f'{prefix}_{set_name}_{idx}']
                
    # Format the data
    rncmpt_list = np.array([batch_metadata_df['Motif_display_ID'].to_list()]).T
    gene_name_list = np.array([batch_metadata_df['Gene Name'].to_list()]).T
    species_list = np.array([batch_metadata_df['Species_display_name'].to_list()]).T
    class_list = np.array([class_list]).T
    prob_list = np.array([prob_list]).T
    kmer_mat = top_kmer_mat.reshape(len(top_kmer_mat) // 3, -1)
    zscore_mat = top_zscore_mat.reshape(len(top_zscore_mat) // 3, -1)
    aligned_kmer_mat = aligned_top_kmer_mat.reshape(len(aligned_top_kmer_mat) // 3, -1)
    summary_mat = np.hstack((rncmpt_list, gene_name_list, species_list, class_list, prob_list,
                             kmer_mat, zscore_mat, aligned_kmer_mat))
    summary_df = pd.DataFrame(summary_mat, index=batch_metadata_df['Hyb_display_ID'], columns=col_list)
    return summary_df


def main(batch):
    """
    Entry point of the program.
    
    """

    # Configure the logging behavior
    logging.set_verbosity(logging.INFO)

    # Get directories
    batch_folder = os.path.join('..', batch)
    rnacompete_metadata_path = 'metadata/rnacompete_metadata.csv'
    probe_metadata_path = 'metadata/probe_metadata.csv'
    probe_intensity_path = os.path.join(batch_folder, 'probe_intensity.csv')

    probe_zscore_path = os.path.join(batch_folder, 'probe_zscore.csv')
    kmer_zscore_path = os.path.join(batch_folder, 'kmer_zscore.csv')
    scatter_folder = os.path.join(batch_folder, 'scatter')
    pwm_folder = os.path.join(batch_folder, 'pwm')
    logo_folder = os.path.join(batch_folder, 'logo')
    feature_path = os.path.join(batch_folder, 'feature.csv')
    summary_path = os.path.join(batch_folder, 'summary.csv')

    # Create directories
    if not os.path.exists(scatter_folder):
        os.makedirs(scatter_folder)
    if not os.path.exists(pwm_folder):
        os.makedirs(pwm_folder)
    if not os.path.exists(logo_folder):
        os.makedirs(logo_folder)

    # Read input files
    rnacompete_metadata_df = pd.read_csv(rnacompete_metadata_path, sep='\t')
    probe_metadata_df = pd.read_csv(probe_metadata_path, sep='\t', index_col=0)
    probe_intensity_df = pd.read_csv(probe_intensity_path, sep='\t', index_col=0)

    # Step 1 | Compute probe Z-scores
    logging.info('STEP 1 | Computing probe Z-scores')
    probe_zscore_df = compute_probe_zscore(probe_intensity_df, probe_metadata_df)
    probe_zscore_df.to_csv(probe_zscore_path, sep='\t')

    # Step 2 | Compute k-mer Z-scores
    logging.info('STEP 2 | Computing k-mer Z-scores')
    kmer_zscore_df = compute_kmer_zscore(probe_zscore_df, probe_metadata_df)
    kmer_zscore_df.to_csv(kmer_zscore_path, sep='\t')

    # Step 3 | Plot scatter plot
    logging.info('STEP 3 | Plotting scatter plots')
    plot_scatter(kmer_zscore_df, scatter_folder)

    # Step 4 | Get top 10 k-mers
    logging.info('STEP 4 | Getting top 10 k-mers')
    top_kmer_mat, top_zscore_mat = get_top_kmer(kmer_zscore_df)

    # Step 5 | Align top 10 k-mers
    logging.info('STEP 5 | Aligning top 10 k-mers')
    aligned_top_kmer_mat = align_top_kmer(top_kmer_mat)

    # Step 6 | Compute PWMs
    logging.info('STEP 6 | Computing PWMs')
    hybid_list = probe_zscore_df.columns.to_list()
    hybid_set_list = kmer_zscore_df.columns.to_list()
    pwm_mat_list = compute_pwm(hybid_set_list, aligned_top_kmer_mat, top_zscore_mat, pwm_folder)

    # Step 7 | Plot logos
    logging.info('STEP 7 | Plotting logos')
    plot_logo(hybid_set_list, pwm_mat_list, logo_folder)

    # Step 8 | Run classifier
    logging.info('STEP 8 | Running classifier')
    feature_df, prob_list = run_classifier(kmer_zscore_df, top_kmer_mat, top_zscore_mat,
                                           pwm_mat_list, pwm_folder)
    class_list = ['Success' if e >= 0.5 else 'Failure' for e in prob_list]
    feature_df.to_csv(feature_path, sep='\t')

    # Step 9 | Get summary
    logging.info('STEP 9 | Get summary')
    batch_metadata_df = rnacompete_metadata_df[rnacompete_metadata_df['Hyb_display_ID'].isin(hybid_list)]
    summary_df = get_summary(batch_metadata_df, top_kmer_mat, top_zscore_mat, aligned_top_kmer_mat,
                             class_list, prob_list)
    summary_df.replace('nan', '').to_csv(summary_path, sep='\t')


if __name__ == '__main__':
    main(sys.argv[1])
