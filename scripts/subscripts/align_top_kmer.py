import numpy as np


def align(kmer_1, kmer_2):
    """
    Align two k-mers.
    
    The two k-mers are shifted left and right against one
    another to find the optimal offset to maximize matching
    nucleotides.
    
    Parameters
    ----------
    kmer_1 : str
        First k-mer.
    kmer_2 : str
        Second k-mer.
        
    Returns
    -------
    max_match : int
        Maximum number of matching nucleotides.
    min_offset : int
        Offset required for maximum matching.

    """
    max_match = 0
    min_offset = 0
    k = len(kmer_1)
    
    # Shift kmer 1 left
    for i in range(k):
        match = 0
        cmp = kmer_1[i:]
        for j in range(len(cmp)):
            if cmp[j] == kmer_2[j]:
                match += 1
        if match > max_match:
            max_match = match
            min_offset = i
    
    # Shift kmer 2 left
    for i in range(k):
        match = 0
        cmp = kmer_2[i:]
        for j in range(len(cmp)):
            if cmp[j] == kmer_1[j]:
                match += 1
        if match > max_match:
            max_match = match
            min_offset = -i
    return max_match, min_offset


def align_top_kmer(top_kmer_mat):
    """
    Perform pairwise alignment for all k-mers against the
    seed k-mer, for all proteins and sets
    
    The seed k-mer is selected as the k-mer with the most
    average matches against the other k-mers, with the least
    average offset.
    
    Parameters
    ----------
    top_kmer_mat : np.ndarray
        Matrix of the top 10 k-mers for the proteins in different sets.
        
    Returns
    -------
    aligned_top_kmer_mat : np.ndarray
        Matrix of the top 10 aligned k-mers for the proteins in different sets.
    
    """
    
    # Iterate over the proteins and sets
    aligned_top_kmer_mat = []
    for kmer_list in top_kmer_mat:

        # Find the seed k-mer
        match_sum_list = []
        offset_sum_list = []

        # Iterate over the k-mers
        for kmer_1 in kmer_list:
            match_sum = 0
            offset_sum = 0

            # Iterate over the k-mers
            for kmer_2 in kmer_list:

                # Perform alignment and get total sum and offset
                match, offset = align(kmer_1, kmer_2)
                match_sum += match
                #offset_sum += abs(offset)
                offset_sum += 0
            match_sum_list.append(match_sum)
            offset_sum_list.append(offset_sum)
        match_sum_list = np.array(match_sum_list)
        offset_sum_list = np.array(offset_sum_list)

        # Select the seed k-mer
        best_idx_list = np.argwhere(match_sum_list == match_sum_list.max()).flatten()
        best_idx = best_idx_list[np.argmin(offset_sum_list[best_idx_list])]

        # Perform alignment using the seed k-mer
        offset_list = []
        for kmer in kmer_list:
            match, offset = align(kmer_list[best_idx], kmer)
            offset_list.append(offset)

        # Format alignments
        aligned_top_kmer_list = []
        min_offset = min(offset_list)
        max_offset = max(offset_list)
        for kmer, offset in zip(kmer_list, offset_list):
            aligned_top_kmer_list.append('-' * (offset - min_offset) + kmer + '-' * (max_offset - offset))
        aligned_top_kmer_mat.append(aligned_top_kmer_list)
    aligned_top_kmer_mat = np.array(aligned_top_kmer_mat)
    return aligned_top_kmer_mat
