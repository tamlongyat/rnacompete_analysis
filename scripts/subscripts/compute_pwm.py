import os

import numpy as np


def save_pwm(hybid_set, pwm_mat, pwm_path):
    """
    Save the PWM as a MEME file.
    
    Parameters
    ----------
    hybid_set : str
        HybID and set.
    pwm_mat : np.ndarray
        PWM.
    pwm_path : str
        Path of the output MEME file. 
    
    """
    
    # Write the header
    line_list = ['MEME version 4\n\n', 'ALPHABET= ACGU\n\n', 'strands: +\n\n',
                 'Background letter frequencies\n', 'A 0.25 C 0.25 G 0.25 T 0.25\n\n',
                 f'MOTIF {hybid_set}\n', 'letter-probability matrix:\n']
    
    # Write the PWM
    for pwm_list in pwm_mat:
        for e in pwm_list:
            line_list += f'{e:06f} '
        line_list += '\n'
    
    # Save the file
    with open(pwm_path, 'w') as f:
        f.writelines(line_list)


def compute_pwm(hybid_set_list, aligned_top_kmer_mat, top_zscore_mat, pwm_folder):
    """
    Compute and save the PWMs for a list of proteins in different sets.
    
    The PWMs only contain positions where no less than half
    of the k-mers are represented.
    
    Parameters
    ----------
    hybid_set_list : list
        List of the HybIDs and their sets.
        For example: HybID00001_a, HybID00001_b, HybID00001_ab.
    aligned_kmer_mat : np.ndarray
        Matrix of the top 10 aligned k-mers for the proteins in different sets.
    top_zscore_mat : np.ndarray
        Matrix of the Z-scores of the top 10 k-mers for the proteins in different sets.
    pwm_folder : str
        Folder to save the PWMs.
    
    Returns
    -------
    pwm_mat_list : list
        List of the PWMs.

    """

    # Iterate over the proteins and sets
    pwm_mat_list = []
    for hybid_set, kmer_list, zscore_list in zip(hybid_set_list,
                                                 aligned_top_kmer_mat,
                                                 top_zscore_mat):

        # Initialize the variables
        mean_zscore = np.mean(zscore_list)
        num_pos = len(kmer_list[0])
        num_kmer = len(kmer_list)
        pfm_mat = np.zeros((num_pos, 4))
        pwm_mat = np.zeros((num_pos, 4))

        # Iterate over the positions
        for pos in range(num_pos):
            
            # Iterate over the k-mers
            for kmer_idx, kmer in enumerate(kmer_list):
                
                # Iterate over the nucleotides
                for nt_idx, nt in enumerate(['A', 'C', 'G', 'U']):
                    if kmer[pos] == nt:
                        pfm_mat[pos, nt_idx] += 1
                        pwm_mat[pos, nt_idx] += zscore_list[kmer_idx]
                    if kmer[pos] == '-':
                        pwm_mat[pos, nt_idx] += mean_zscore / 4

        # Trim the PWM
        pwm_mat = pwm_mat[pfm_mat.sum(axis=1) / 10 >= 0.5]

        # Add pseudo-count and convert to fraction
        pwm_mat += 1
        pwm_mat /= pwm_mat.sum(axis=1)[:, None]
        
        # Save the PWM
        pwm_mat_list.append(pwm_mat)
        pwm_path = os.path.join(pwm_folder, f'{hybid_set}.txt')
        save_pwm(hybid_set, pwm_mat, pwm_path)
    return pwm_mat_list
