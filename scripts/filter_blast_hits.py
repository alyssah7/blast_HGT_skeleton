"""
Created on Tue Sep  2 17:14:41 2025

@author: Alyssa Mitchell
"""
import numpy as np
import pandas as pd
import os
import sys
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='''\
                            filters the all-by-all blast table for length and contig quality
                               ''',
                               epilog="Questions or comments? --> Alyssah@mit.edu")
parser.add_argument("-s", dest="summary_table", type=str, help="Path to the summary table (.tsv) for all BLAST runs",required=True,action='store')
parser.add_argument("-l", dest="length", type=int, help="minimum length of BLAST hit to consider",required=True)
parser.add_argument("-c", dest="cov", type=int, help="minimum SPAdes kmer coverage to consider contig",required=True)
parser.add_argument("-o", dest="outpath", type=str, help="Path to the output/filtered summary table (.tsv)",required=True)
parser.add_argument("-g", dest="contig_info", type=str, help="path to csv file with each genome's longest contig info",required=True)
parser.add_argument("-e", dest="lenthresh", type=int, help="Threshold for minimum length of longest contig to consider HGT in this isolate",required=True)
parser.add_argument("-k", dest="covthresh", type=int, help="Threshold for minimum kmer coverage of longest contig to consider HGT in this isolate",required=True)
parser.add_argument("-f", dest="fraccov", type=float, help="Threshold for minimum fraction of largest contig kmer coverage to consider HGT on this contig",required=True)

args = parser.parse_args()

def execute_filter(in_path, len_thresh, cov_thresh, out_path, genome_contig, l_thresh, k_thresh, frac_thresh):
    # test
    # df = pd.read_csv('Desktop/blast_AHM_v0007_A01__AHM_v0007_A03_out.tsv', sep='\t', header=None).rename(columns={0: 'query', 1: 'subject', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore'})
    # df = pd.read_csv('Desktop/blast_AHM_v0007_B09__AHM_v0009_B10_out.tsv', sep='\t', header=None).rename(columns={0: 'query', 1: 'subject', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore'})

    df = pd.read_csv(in_path, sep='\t', header=None).rename(columns={0: 'q_contig', 1: 's_contig', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore', 12: 'query', 13: 'subject'})

    # only consider matches with sufficient length (input as arg)
    df = df[df['length'] >= len_thresh]
    
    # only keep hits in contigs with good kmer coverage from SPAdes
    df['qcov'] = [float(q.split('_cov_')[1]) for q in df['q_contig']]
    df['scov'] = [float(q.split('_cov_')[1]) for q in df['s_contig']]
    df = df[(df['qcov'] >= cov_thresh) & (df['scov'] >= cov_thresh)]

    assembly_df = pd.read_csv(genome_contig)

    # test
    # assembly_df = pd.read_csv('Desktop/longest_contigs_all.csv')
    # df = pd.read_csv('Desktop/filtered_blast_AHM_v0006_F09.tsv', sep=',', header=None).rename(columns={0: 'q_contig', 1: 's_contig', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore', 12: 'query', 13: 'subject', 14: 'qcov', 15: 'scov'})
    df['query'] = [q.split('__')[0].split('blast_')[1] for q in df['query']]
    df['subject'] = [s.split('__')[1].split('_out')[0] for s in df['subject']]
    
    # filter out genomes with poor assemblies (by length of longest contig and kmer cov of that contig)
    filt_assemblies = assembly_df[assembly_df['longest_length'] >= l_thresh]
    filt_assemblies = filt_assemblies[filt_assemblies['longest_cov'] >= k_thresh]
    
    filt_df = df[df['query'].isin(filt_assemblies['sample']) & df['subject'].isin(filt_assemblies['sample'])]
    
    # now filter out any hits on contigs less than n% coverage of largest contig for that genome
    filt_df['qcov_frac'] = filt_df['qcov']/[filt_assemblies[filt_assemblies['sample'] == q]['longest_cov'].values[0] for q in filt_df['query']]
    filt_df['scov_frac'] = filt_df['scov']/[filt_assemblies[filt_assemblies['sample'] == q]['longest_cov'].values[0] for q in filt_df['subject']]
    
    filt_df = filt_df[(filt_df['qcov_frac'] >= frac_thresh) & (filt_df['scov_frac'] >= frac_thresh)]
    
    filt_df.to_csv(out_path, header=False, index=False)

    return

# %%
if __name__ == "__main__":
    execute_filter(args.summary_table, args.length, args.cov, args.outpath, args.contig_info, args.lenthresh, args.covthresh, args.fraccov)
