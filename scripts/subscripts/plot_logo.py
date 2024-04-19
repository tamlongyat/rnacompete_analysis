import os

import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_logo(hybid_set_list, pwm_mat_list, logo_folder):
    """
    Plot and save the logos for a list of PWMs.
    
    Parameters
    ----------
    hybid_set_list : list
        List of the HybIDs and their sets.
        For example: HybID00001_a, HybID00001_b, HybID00001_ab.
    pwm_mat_list : list
        List of the PWMs.
    logo_folder : str
        Folder to save the logos.

    """
    
    # Define the color scheme
    color_scheme = {
        'A': '#00CC00', 
        'C': '#0000CC',
        'G': '#FFB302',
        'U': '#CC0001',
    }
    
    # Iterate over the proteins and sets
    for hybid_set, pwm_mat in zip(hybid_set_list, pwm_mat_list):
        
        # Convert PWM to bits
        bit_list = np.clip(2 + np.sum(pwm_mat * np.log2(pwm_mat), axis=1), 0, None)
        logo_mat = pwm_mat * bit_list[:, None]
        
        # Plot logo
        logomaker.Logo(pd.DataFrame(logo_mat, columns=['A', 'C', 'G', 'U']),
                       color_scheme=color_scheme,
                       show_spines=False)
        plt.xticks([])
        plt.yticks([])
        
        # Save logo
        logo_path = os.path.join(logo_folder, f'{hybid_set}.png')
        plt.savefig(logo_path)
        plt.close()
