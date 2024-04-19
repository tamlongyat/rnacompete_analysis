import argparse
import glob
import os
import shutil
import sys

from absl import logging
import numpy as np
import pandas as pd


def generate_subpage(hybid, summary_list):
    """
    Generate the subpage HTML file.
    
    Parameters
    ----------
    protein_id : str
        Protein ID
    summary_list : list
        Row of the summary dataframe that corresponds to the protein.
    
    """
    rncmpt = summary_list['rncmpt']
    gene_name = summary_list['gene_name']
    species = summary_list['species']
    cls = summary_list['class']
    prob = summary_list['prob']
    
    # Write header
    html_str = f'<h1>{hybid} | {rncmpt} | {gene_name} | {species} | {cls} ({prob:.2f})</h1>\n'
    
    # Write scatter plot
    html_str += '<h2>7-mer Scatter plot</h2>\n'
    html_str += '<img src="scatter.png">\n'
    html_str += '<hr>\n'
    
    # Write logos
    logo_df = pd.DataFrame({'Set A': ['<img src="logo_a.png" height=80>'],
                            'Set B': ['<img src="logo_b.png" height=80>'],
                            'Set AB': ['<img src="logo_ab.png" height=80>']})
    html_str += '<h2>Motifs</h2>\n'
    html_str += logo_df.to_html(index=False, justify='center', escape=False)
    html_str += '<hr>\n'
    
    # Write K-mers and Z-scores
    summary_mat = summary_list[5:95].to_numpy().reshape(-1, 10)
    col_list = ['Set A K-mer', 'Set A Z-score', 'Set A Aligned', '⬛⬛',
                'Set B K-mer', 'Set B Z-score', 'Set B Aligned', '⬛⬛',
                'Set AB K-mer', 'Set AB Z-score', 'Set AB Aligned', '⬛⬛',]
    summary_mat = np.array([[f'<tt>{e}</tt>' for e in f] for f in summary_mat])
    kmer_df = pd.DataFrame(np.vstack((summary_mat[[0, 3, 6]], ['⬛⬛'] * 10,
                                      summary_mat[[1, 4, 7]], ['⬛⬛'] * 10,
                                      summary_mat[[2, 5, 8]], ['⬛⬛'] * 10)).T,
                           columns=col_list)
    html_str += '<h2>Top 7-mers and Z-scores</h2>\n'
    html_str += kmer_df.to_html(index=False, justify='center', escape=False)
    html_str += '<hr>\n'
    html_str += '<a href="../index.html">Go back</a>'
    with open(f'../html/{hybid}/{hybid}.html', 'w') as f:
        f.write(html_str)
        

def generate_index(summary_df):
    """
    Generate the index HTML file.
    
    Parameters
    ----------
    summary_df : pd.DataFrame
        Dataframe containing the summary of the selected proteins.
    
    """
    index_df = pd.DataFrame()
    index_df['Hyb ID'] = [f'<a href="{e}/{e}.html">{e}</a>' for e in summary_df.index]
    index_df['RNAcompete ID'] = summary_df['rncmpt'].to_list()
    index_df['Gene name'] = summary_df['gene_name'].to_list()
    index_df['Species'] = summary_df['species'].to_list()
    classifier_list = []
    for cls, prob in zip(summary_df['class'], summary_df['prob']):
        if cls == 'Failure':
            classifier_list.append(f'<span style="color:red">{cls} ({prob:.2f})</span>')
        else:
            classifier_list.append(f'<span style="color:green">{cls} ({prob:.2f})</span>')
    index_df['Classification'] = classifier_list
    index_df['Logo'] = [f'<img src="{e}/logo_ab.png" height=80>' for e in summary_df.index]
    
    # Convert table to HTML
    html_str = index_df.to_html(index=False, na_rep='', justify='center', escape=False)
    with open('../html/index.html', 'w') as f:
        f.write(html_str)


def main(selection, nohtml):
    """
    Entry point of the program.
    
    """
    
    # Configure the logging behavior
    logging.set_verbosity(logging.INFO)
    
    # Create output folder
    if not nohtml and not os.path.exists('../html'):
        os.makedirs('../html')

    # Read selected proteins
    if selection is not None:
        with open(selection) as f:
            selected_hybid_list = f.read().splitlines()
        logging.info(f'{len(selected_hybid_list)} proteins selected')
    
    # Iterate over batches
    summary_df_list = []
    batch_folder_list = sorted(glob.glob('../HybID*'))
    for batch_folder in batch_folder_list:
        tmp_list = glob.glob(os.path.join(batch_folder, 'logo', '*ab.png'))
        hybid_list = sorted([os.path.basename(e).split('_')[0] for e in tmp_list])
        
        # Get selected proteins
        if selection is not None:
            hybid_list = sorted(list(set(hybid_list).intersection(set(selected_hybid_list))))

        # Read summary
        if len(hybid_list) > 0:
            summary_df = pd.read_csv(os.path.join(batch_folder, 'summary.csv'), sep='\t', index_col=0)
            logging.info(f'Generating HTML files for proteins in {batch_folder}')

        # Iterate over proteins
        for hybid in hybid_list:
        
            # Read summary
            summary_list = summary_df.loc[hybid]
            summary_df_list.append(summary_list)

            if not nohtml:
            
                # Create folder
                if not os.path.exists(f'../html/{hybid}'):
                    os.makedirs(f'../html/{hybid}')

                # Copy logos and scatter plots
                shutil.copyfile(os.path.join(batch_folder, 'logo', f'{hybid}_a.png'),
                                f'../html/{hybid}/logo_a.png')
                shutil.copyfile(os.path.join(batch_folder, 'logo', f'{hybid}_b.png'),
                                f'../html/{hybid}/logo_b.png')
                shutil.copyfile(os.path.join(batch_folder, 'logo', f'{hybid}_ab.png'),
                                f'../html/{hybid}/logo_ab.png')

                shutil.copyfile(os.path.join(batch_folder, 'scatter', f'{hybid}.png'),
                                f'../html/{hybid}/scatter.png')

                # Generate HTML file
                generate_subpage(hybid, summary_list)

    # Generate index file
    summary_df = pd.DataFrame(summary_df_list)
    if not nohtml:
        logging.info(f'Generating index HTML file')
        generate_index(summary_df)
    
    # Save summary
    summary_df.to_csv('../summary.csv', sep='\t')


if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate HTML reports for all experiments.')
    parser.add_argument('selection', type=str, nargs='?', default=None,
                        help='Filename of a file containing the names of selected experiments.')
    parser.add_argument('--nohtml', action='store_true',
                        help='Whether not to generate HTML reports.')
    args = parser.parse_args()
    main(args.selection, args.nohtml)
