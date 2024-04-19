import os
import subprocess

import numpy as np
import pandas as pd
import pickle
from scipy import stats


def run_tomtom(pwm_folder, hybid, set_1, set_2):
    """
    Compute the similarity between two motifs.
    
    Parameters
    ----------
    pwm_folder : str
        Folder containing the PWMs.
    hybid : str
        HybID
    set_1 : str
        Set A or set B
    set_2 : str
        Set A or set B
    
    Returns
    -------
    sim
        Negative log10(p-value).
    
    """
            
    # Set paths
    pwm_1_path = os.path.join(pwm_folder, f'{hybid}_{set_1}.txt')
    pwm_2_path = os.path.join(pwm_folder, f'{hybid}_{set_2}.txt')
            
    # Run TOMTOM
    cmd_1_list = ['tomtom', '-min-overlap', '4', '-norc', '-text', '-thresh', '1', pwm_1_path, pwm_2_path]
    cmd_2_list = ['sed', '-n', '2p']
    ps = subprocess.Popen(cmd_1_list, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    output = subprocess.check_output(cmd_2_list, stdin=ps.stdout).decode()
    sim = -np.log10(float(str(output).split('\t')[3]))
    return sim
  

def generate_feature(kmer_zscore_df, top_kmer_mat, top_zscore_mat, pwm_mat_list,
                     pwm_folder):
    """
    Generate features for all proteins for the classifier.
    
    Parameters
    ----------
    kmer_zscore_df : pd.DataFrame
        DataFrame containing the k-mer Z-scores for different proteins and sets.
    top_kmer_mat : np.ndarray
        Matrix of the top 10 k-mers for the proteins in different sets.
    top_zscore_mat : np.ndarray
        Matrix of the Z-scores of the top 10 k-mers for the proteins in different sets.
    pwm_mat_list : list
        List of the PWMs.
    pwm_folder : str
        Folder containing the PWMs.
    
    Returns
    -------
    feature_df : np.ndarray
        Feature DataFrame.
    
    """
    
    # Define artifacts and column names
    artifact_list = ['AAAUAAA', 'AGACCC', 'AGACGGG', 'CCCC', 'CGGAGG',
                     'GACUCAC', 'GACUGCC', 'GAGUC', 'GGGAGC', 'GGGG',
                     'GGGGAGC', 'GGGGCC', 'GGGCC', 'GGGGGGG', 'AGGCCC',
                     'AGGGCC', 'UUUUUUU', 'AGGGC', 'CCCCC', 'CCCCCC',
                     'GGGGG', 'GGGGGG', 'CCCCCCC', 'GGCGG', 'AGAGCC',
                     'GUCGU']
    
    # Initialize empty result
    feature_mat = []
    
    # Iterate over the protein IDs
    hybid_list = np.unique([e.split('_')[0] for e in kmer_zscore_df.columns])
    for idx, hybid in enumerate(hybid_list):

        # Read the Z-scores
        zscore_a_list = kmer_zscore_df[f'{hybid}_a'].dropna()
        zscore_b_list = kmer_zscore_df[f'{hybid}_b'].dropna()
        zscore_ab_list = kmer_zscore_df[f'{hybid}_ab'].dropna()
        
        # Read the top k-mers
        top_kmer_a_list = top_kmer_mat[3 * idx]
        top_kmer_b_list = top_kmer_mat[3 * idx + 1]
        top_kmer_ab_list = top_kmer_mat[3 * idx + 2]
        top_kmer_list = list(top_kmer_a_list) + list(top_kmer_b_list)
        
        # Read the top Z-scores
        top_zscore_a_list = top_zscore_mat[3 * idx]
        top_zscore_b_list = top_zscore_mat[3 * idx + 1]
        top_zscore_ab_list = top_zscore_mat[3 * idx + 2]
        
        # Read the PWMs
        pwm_a_mat = pwm_mat_list[3 * idx]
        pwm_b_mat = pwm_mat_list[3 * idx + 1]
        pwm_ab_mat = pwm_mat_list[3 * idx + 2]

        # Compute the Z-score statistics
        pearson = stats.pearsonr(zscore_a_list, zscore_b_list)[0]
        num_intersect = len(set(top_kmer_a_list).intersection(set(top_kmer_b_list)))
        skew_a = stats.skew(zscore_a_list)
        skew_b = stats.skew(zscore_b_list)
        kurt_a = stats.kurtosis(zscore_a_list, fisher=False)
        kurt_b = stats.kurtosis(zscore_b_list, fisher=False)
        artifact_count_list = [np.sum([artifact in e for e in top_kmer_list]) for artifact in artifact_list]
        artifact_count_sum = sum(artifact_count_list)

        # Compute the information content
        ic_a = 14 + (pwm_a_mat * np.log2(pwm_a_mat)).sum()
        ic_b = 14 + (pwm_b_mat * np.log2(pwm_b_mat)).sum()
        
        # Compute the motif similarity
        motif_sim_ab = run_tomtom(pwm_folder, hybid, 'a', 'b')
        motif_sim_ba = run_tomtom(pwm_folder, hybid, 'b', 'a')

        # Combine features
        feature_list = [pearson, num_intersect, *top_zscore_a_list, *top_zscore_b_list,
                        skew_a, skew_b, kurt_a, kurt_b,
                        *artifact_count_list, artifact_count_sum,
                        ic_a, ic_b, motif_sim_ab, motif_sim_ba,
                        top_zscore_ab_list[0]]
        feature_mat.append(feature_list)
    
    # Convert matrix into DataFrame
    feature_mat = np.array(feature_mat)
    col_list = ['zscoreCorrelation', 'top10overlap',
                'aZ1', 'aZ2', 'aZ3', 'aZ4', 'aZ5', 'aZ6', 'aZ7', 'aZ8', 'aZ9', 'aZ10',
                'bZ1', 'bZ2', 'bZ3', 'bZ4', 'bZ5', 'bZ6', 'bZ7', 'bZ8', 'bZ9', 'bZ10',
                'aSkewness', 'bSkewness', 'aKurtosis', 'bKurtosis',
                *artifact_list, 'total_artifacts',
                'aIC', 'bIC', 'AtoB', 'BtoA', 'setAB.zscore']
    feature_df = pd.DataFrame(feature_mat, index=hybid_list, columns=col_list)
    return feature_df


def run_classifier(kmer_zscore_df, top_kmer_mat, top_zscore_mat, pwm_mat_list, pwm_folder):
    """
    Entry point of the program
    
    """

    # Get features
    feature_df = generate_feature(kmer_zscore_df, top_kmer_mat, top_zscore_mat,
                                  pwm_mat_list, pwm_folder)
    
    # Load model
    with open('classifier/SS_model.sav', 'rb') as f:
        ss = pickle.load(f)
    with open('classifier/LR_model.sav', 'rb') as f:
        lr = pickle.load(f)
    
    # Run model
    feature_mat = ss.transform(feature_df.to_numpy())
    prob_list = lr.predict_proba(feature_mat)[:, 1]
    return feature_df, prob_list
