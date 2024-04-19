import numpy as np
from scipy.stats import median_abs_deviation


def compute_probe_zscore(probe_intensity_df, probe_metadata_df):
    """Compute the probe Z-scores from a probe intensity DataFrame.
    
    The Z-score computation is a three-step process:
    - Step 1: Probe intensities are normalized by the median and IQR column-wise. They are then rescaled to
    the geometric mean of the medians and IQRs.
    - Step 2: Probe intensities are normalized by the median and MAD row-wise.
    - Step 3: Probe intensities are normalized by the median and MAD column-wise.
    
    Parameters
    ----------
    probe_intensity_df : DataFrame
        Probe intensity DataFrame.
    probe_metadata_df : DataFrame
        Metadata DataFrame which contains the sequences and sets of probes.
        
    Returns
    -------
    step_3_df : DataFrame
        Probe Z-score DataFrame.
    
    """
    
    # Subset probes in the probe metadata DataFrame
    probe_id_list = probe_metadata_df.index
    probe_intensity_df = probe_intensity_df.loc[probe_id_list]
    
    # Set flagged values to NaN
    med_df = probe_intensity_df.iloc[:, ::2]
    flag_mat = probe_intensity_df.iloc[:, 1::2].to_numpy()
    step_0_df = med_df.copy()
    step_0_df[flag_mat > 0] = np.NaN
    
    # Step 1
    median_1_list = np.nanmedian(step_0_df, axis=0)
    iqr_1_list = (np.nanquantile(step_0_df, 0.75, axis=0, method='hazen')
                  - np.nanquantile(step_0_df, 0.25, axis=0, method='hazen'))
    step_1_df = (step_0_df - median_1_list) / iqr_1_list
    new_median_1 = np.exp(np.nanmean(np.log(median_1_list)))
    new_iqr_1 = np.exp(np.nanmean(np.log(iqr_1_list)))
    step_1_df = step_1_df * new_iqr_1 + new_median_1
    
    # Step 2
    median_2_list = np.expand_dims(np.nanmedian(step_1_df, axis=1), axis=1)
    mad_2_list = np.expand_dims(median_abs_deviation(step_1_df, axis=1, nan_policy='omit'), axis=1)
    step_2_df = (step_1_df - median_2_list) / (1.4826 * mad_2_list)
    
    # Step 3
    median_3_list = np.nanmedian(step_2_df, axis=0)
    mad_3_list = median_abs_deviation(step_2_df, axis=0, nan_policy='omit')
    probe_zscore_df = (step_2_df - median_3_list) / (1.4826 * mad_3_list)
    return probe_zscore_df
