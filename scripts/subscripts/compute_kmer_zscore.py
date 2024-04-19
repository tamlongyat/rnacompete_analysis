import itertools

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, find


def compute_kmer_count(probe_metadata_df, k=7):
    """Compute the k-mer count per probe in the probe metadata.
    
    Parameters
    ----------
    probe_metadata_df : 
        Metadata DataFrame which contains the sequences and sets of probes.
    k : int
        Length of the k-mer,
    
    Returns
    -------
    probe_id_list : list
        List of all probe IDs.
    probe_set_list : np.ndarray
        List of all probe sets
    kmer_list : list
        List of all k-mers.
    kmer_count_mat : Sparse CSC matrix (num_probe, num_kmer)
        Matrix indicating the k-mer count per probe.
    
    """
    
    # Get the list of probe IDs, probe sequences and sets
    probe_id_list = probe_metadata_df.index
    probe_seq_list = probe_metadata_df['seq']
    probe_set_list = probe_metadata_df['set'].to_numpy()
    
    # Generate the list of kmers
    nuc_list = ['A', 'C', 'G', 'U']
    kmer_list = [''.join(i) for i in itertools.product(nuc_list, repeat=k)]
    kmer_dict = {k:v for v, k in enumerate(kmer_list)}
    
    # Initialize the k-mer count matrix
    kmer_count_mat = np.zeros((len(probe_id_list), len(kmer_list)), dtype=np.int8)
    
    # Fill in the k-mer count matrix
    for probe_seq_idx, probe_seq in enumerate(probe_seq_list):
        for start_idx in range(len(probe_seq) - k + 1):
            kmer_idx = kmer_dict[probe_seq[start_idx:start_idx + k]]
            kmer_count_mat[probe_seq_idx, kmer_idx] += 1
    kmer_count_mat = csc_matrix(kmer_count_mat)    
    
    return probe_id_list, probe_set_list, kmer_list, kmer_count_mat


def compute_kmer_zscore_set_kmer(probe_zscore_kmer_mat, trim_ratio=0.025):
    """Compute the k-mer Z-scores from a subset probe Z-score matrix for a given k-mer.
    
    For each RNCMPT IDs, the k-mer Z-score is computed as the trimmed mean of all probe Z-scores
    excluding NaNs, with the top 2.5% and bottom 2.5% dropped.
    
    Parameters
    ----------
    probe_zscore_kmer_mat : np.ndarray
        Probe Z-score Matrix for probes containing the given k-mer.
    trim ratio: float
        Ratio for the top and bottom Z-scores to be trimmed.
    
    Returns
    -------
    kmer_zscore_kmer_list : np.ndarray
        List of Z-scores for all RNCMPT IDs for the given k-mer.
    
    """
    
    # Trimmed mean function (imitates MATLAB rounding behavior)
    def trim_mean_round_off(a, proportiontocut):
        nobs = len(a)
        k = proportiontocut * nobs
        lowercut = int(np.ceil(k - 0.5))
        uppercut = nobs - lowercut
        atmp = np.partition(a, (lowercut, uppercut - 1))
        return np.mean(atmp[lowercut:uppercut])
        
    # Initialize Z-score list for the kmer
    num_rncmpt_id = probe_zscore_kmer_mat.shape[1]
    kmer_zscore_kmer_list = np.zeros(num_rncmpt_id)
        
    # Iterate over RNCMPT IDs
    for rncmpt_id_idx in range(num_rncmpt_id):
            
        # Step 1
        probe_zscore_list = probe_zscore_kmer_mat[:, rncmpt_id_idx]
            
        # Step 2
        probe_zscore_list = probe_zscore_list[~np.isnan(probe_zscore_list)]
            
        # Step 3
        if len(probe_zscore_list) == 0:
            kmer_zscore_kmer_list[rncmpt_id_idx] = np.NaN
        else:
            kmer_zscore_kmer_list[rncmpt_id_idx] = trim_mean_round_off(probe_zscore_list, trim_ratio)
        
    return kmer_zscore_kmer_list


def compute_kmer_zscore_set(probe_zscore_mat, kmer_count_mat):
    """Compute the k-mer Z-scores from a probe Z-score matrix for all k-mers.
    
    For each k-mer, the probe Z-score matrix is subset for probes which contain the k-mer.
    Each of these subset probe Z-score matrices is then processed in parallel to compute
    the k-mer Z-score for all RNCMPT IDs.
    
    The resulting k-mere Z-scores are then normalized by the mean and standard deviation
    column-wise.
    
    Parameters
    ----------
    probe_zscore_mat : np.ndarray
        Probe Z-score matrix.
    kmer_count_mat : Sparse CSC matrix (num_probe, num_kmer)
        Matrix indicating the k-mer count per probe.
    
    Returns
    -------
    kmer_zscore_mat : np.ndarray
        k-mer Z-score Matrix
    
    """
    
        
    # Set probe Z-scores below zero to zero
    probe_zscore_mat[probe_zscore_mat < 0] = 0
    
    # Initialize the k-mer Z-score matrices
    num_kmer = kmer_count_mat.shape[1]
    num_rncmpt_id = probe_zscore_mat.shape[1]
    kmer_zscore_mat = np.zeros((num_kmer, num_rncmpt_id))
    
    # Get list of matrices for probes containing different kmers
    probe_zscore_kmer_mat_list = []
    for kmer_idx in range(num_kmer):
        probe_idx_list = find(kmer_count_mat[:, kmer_idx])[0]
        probe_zscore_kmer_mat_list.append(probe_zscore_mat[probe_idx_list])
    
    # Compute the k-mer Z-scores by parallel processing
    kmer_zscore_kmer_list_list = [None] * num_kmer
    for kmer_idx in range(num_kmer):
        kmer_zscore_kmer_list_list[kmer_idx] = compute_kmer_zscore_set_kmer(probe_zscore_kmer_mat_list[kmer_idx])
    kmer_zscore_mat = np.array(kmer_zscore_kmer_list_list)
    
    # Normalize the data by mean and standard deviation column-wise
    kmer_zscore_mat = (kmer_zscore_mat - np.nanmean(kmer_zscore_mat, axis=0)) / np.nanstd(kmer_zscore_mat, axis=0)
    
    return kmer_zscore_mat


def compute_kmer_zscore(probe_zscore_df, probe_metadata_df):
    """Compute the k-mer Z-scores from a probe Z-score DataFrame for set A, set B and set AB.
    
    Parameters
    ----------
    probe_zscore_df : DataFrame
        Probe Z-score DataFrame.
    probe_metadata_df : DataFrame
        Metadata DataFrame which contains the sequences and sets of probes.
    
    Returns
    -------
    kmer_zscore_a_df : DataFrame
        k-mer Z-score DataFrame for probes in set A.
    kmer_zscore_b_df : DataFrame
        k-mer Z-score DataFrame for probes in set B.
    kmer_zscore_ab_df : DataFrame
        k-mer Z-score DataFrame for probes in set AB.
    
    """
    
    # Get k-mer counts from probe sequences
    probe_id_list, probe_set_list, kmer_list, kmer_count_mat = compute_kmer_count(probe_metadata_df)

    # Subset probes in the probe metadata DataFrame
    probe_zscore_df = probe_zscore_df.loc[probe_id_list]
    
    # Split the probe Z-score DataFrame into set A and B
    probe_zscore_a_mat = probe_zscore_df[probe_set_list == 'SetA'].to_numpy()
    probe_zscore_b_mat = probe_zscore_df[probe_set_list == 'SetB'].to_numpy()
    probe_zscore_ab_mat = probe_zscore_df.to_numpy()
    
    # Split the k-mer count matrix into set A and B
    kmer_count_a_mat = kmer_count_mat[probe_set_list == 'SetA']
    kmer_count_b_mat = kmer_count_mat[probe_set_list == 'SetB']
    kmer_count_ab_mat = kmer_count_mat
    
    # Compute the k-mer Z-score for set A, set B and both
    kmer_zscore_a_mat = compute_kmer_zscore_set(probe_zscore_a_mat, kmer_count_a_mat)
    kmer_zscore_b_mat = compute_kmer_zscore_set(probe_zscore_b_mat, kmer_count_b_mat)
    kmer_zscore_ab_mat = compute_kmer_zscore_set(probe_zscore_ab_mat, kmer_count_ab_mat)
    
    # Combine the matrices
    hybid_list = probe_zscore_df.columns
    kmer_zscore_mat = np.empty((len(kmer_list), 3 * len(hybid_list)))
    col_list = []
    for i, hybid in enumerate(hybid_list):
        kmer_zscore_mat[:, 3 * i] = kmer_zscore_a_mat[:, i]
        kmer_zscore_mat[:, 3 * i + 1] = kmer_zscore_b_mat[:, i]
        kmer_zscore_mat[:, 3 * i + 2] = kmer_zscore_ab_mat[:, i]
        col_list += [f'{hybid}_a', f'{hybid}_b', f'{hybid}_ab']
    
    # Convert the matrix to dataframe
    kmer_zscore_df = pd.DataFrame(kmer_zscore_mat, index=kmer_list, columns=col_list)
    return kmer_zscore_df
