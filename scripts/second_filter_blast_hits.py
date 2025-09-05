"""
Created on Th Sep 4 2025

@author: Alyssa Mitchell
"""
import numpy as np
import pandas as pd
# import os
# import sys
import argparse
import matplotlib.pyplot as plt
from datetime import datetime

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='''\
                            filters the all-by-all blast table for length and contig quality
                               ''',
                               epilog="Questions or comments? --> Alyssah@mit.edu")
parser.add_argument("-i", dest="filt_blast", type=str, help="path to tsv file with filtered blast hits, concatenated this time",required=True)
parser.add_argument("-o", dest="outpath", type=str, help="Path to dir to write outputs to",required=True)
parser.add_argument("-m", dest="meta", type=str, help="Path to metadata .csv file for isolates (to connect to source subject)",required=True)
parser.add_argument("-l", dest="eventsthresh", type=int, help="Threshold for maximum shared genome events across a pair of different-strain isolates",required=True)

args = parser.parse_args()


# def second_filters(in_path, genome_contig, out_dir, l_thresh, k_thresh, frac_thresh, metadata):
def second_filters(in_path, out_dir, metadata, shared_events_thresh):

    print(f'{datetime.now()} : importing concatenated, filtered data')
    filt_df = pd.read_csv(in_path, sep=',', header=None).rename(columns={0: 'q_contig', 1: 's_contig', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore', 12: 'query', 13: 'subject', 14: 'qcov', 15: 'scov', 16: 'qcov_frac', 17: 'scov_frac'})

    print(f'{datetime.now()} : parsing isolate metadata')
    # import isolate metadata and add column to filt_df for within-subject comparison (binary)
    metadata_df = pd.read_csv(metadata, sep=',')
    #test
    # metadata_df = pd.read_csv('Desktop/2025-04_combined_isolate_metadata.csv', sep=',')
    filt_df['q_source'] = [metadata_df[metadata_df['isolate_ID'] == q]['subject'].values[0] for q in filt_df['query']]
    filt_df['s_source'] = [metadata_df[metadata_df['isolate_ID'] == s]['subject'].values[0] for s in filt_df['subject']]
    filt_df['within-person'] = [q == s for q,s in zip(filt_df['q_source'], filt_df['s_source'])]
    
    print(f'{datetime.now()} : counting HGT events')
    # count # of HGT events to separate and remove within-strain events
    event_counts_pairwise = filt_df.groupby(['query', 'subject']).agg(count=('query', 'size'), total_length=('length', 'sum'), within_person=('within-person', 'first')).reset_index()
    
    '''
    print(f'{datetime.now()} : plotting')
    fig, ax = plt.subplots()
    zoomed_total_events = event_counts_pairwise[event_counts_pairwise['count'] < 50]['count']
    # ax.hist(event_counts_pairwise['count'], bins=100, color='darkmagenta')
    ax.hist(zoomed_total_events, bins=50, color='darkmagenta')
    plt.xlabel('# of HGT events in pair')
    plt.ylabel('# of isolate pairs (excluding pairs without any HGT)')
    fig.savefig(out_dir+'/hist_num_HGT_events_by_pair_zoom.png', dpi=300)
    
    
    fig, ax = plt.subplots()
    zoomed_total_lengths = event_counts_pairwise[event_counts_pairwise['total_length'] < 1000000]['total_length']
    # ax.hist(event_counts_pairwise['total_length'], bins=100, color='darkmagenta')
    ax.hist(zoomed_total_lengths, bins=100, color='darkmagenta')
    plt.xlabel('total length of all HGT events in pair')
    plt.ylabel('# of isolate pairs (excluding pairs without any HGT)')
    fig.savefig(out_dir+'/hist_len_HGT_events_by_pair_zoom.png', dpi=300)
    '''
    
    # Merge back into original dataframe
    filt_df = filt_df.merge(
        event_counts_pairwise[['query', 'subject', 'count', 'total_length']],
        on=['query', 'subject'],
        how='left'
    )
    
    print(f'Number of isolate pairs BEFORE filtering for same strain: {len(filt_df)}')
    # Filter to remove suspected same-strain pairs (where they share >10% genome at this high pident)
    filt_df = filt_df[filt_df['count'] <= shared_events_thresh]
    print(f'Number of isolate pairs AFTER filtering for same strain: {len(filt_df)}')

    print(f'{datetime.now()} : saving data structures')
    # save final blast table after filtering
    filt_df.to_csv(out_dir+'/RE_filtered_blast_all.tsv', index=False)
    
    # save isolates to go through genomad
    isolates_of_interest = list(filt_df['query'])
    isolates_of_interest.extend(list(filt_df['subject']))
    isolates_of_interest = np.unique(isolates_of_interest)
    with open(out_dir+'/isolates_of_interest.txt', 'w') as file:
        for iso in isolates_of_interest:
            file.write(str(iso) + '\n')
    
        
    return

#%%
if __name__ == "__main__":
    second_filters(args.filt_blast, args.outpath, args.meta, args.eventsthresh)

