import pandas as pd
import os
import re
from Bio.SeqRecord import SeqRecord
import numpy as np
import pybedtools
import time
import matplotlib.pyplot as plt
import sys
import subprocess
import shutil
import glob
import scipy.stats as stats
import statsmodels as sms
import statsmodels.sandbox.stats.multicomp
import matplotlib
from sklearn.externals.joblib import Parallel, delayed
import itertools as it
from scipy.signal import argrelextrema
import scipy


def run_command(command):
    print('\nRunnning now!\n')
    print(command)
    output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    print('\nDone\nWith ouput:\n%s' % output)
    print(output)
    
def hello_print(command_list):
    print(command_list)

def faidx_genome(genome_fn, contig_fn):
    '''Generates a genome file used for samtools and such if not already present.
    Input: 
    Abspath filename of output genome_file.
    Abspath filename for genome contig file.'''
    if not os.path.exists(genome_fn):
        sammtools_command = 'samtools faidx %s' % contig_fn
        run_command(samtools_command)
        genome_file_command = 'cat %s.fai | sort -k1,1n | cut -f 1,2 > %s' %(contig_fn, genome_fn)
        run_command(genome_file_command)
    else:
        print('%s already exists!' % genome_fn)

def sort_gff_file(gff_fn):
    '''Sorts gff file and returns the Abspath to the sorted gff file.'''
    gff3sport_pl = '/home/benjamin/genome_assembly/PST79/FALCON/p_assemblies/v9_1/Pst_104E_v12/get_homologues/gff3sort/gff3sort.pl'
    if gff_fn.endswith('gff'):
        out_fn = gff_fn.replace('.gff', '.sorted.gff')
    elif gff_fn.endswith('gff3'):
        out_fn = gff_fn.replace('.gff3', '.sorted.gff3')
    else:
        print('Please make sure your gff %s file ends with gff or gff3' % gff_fn.split['/'][-1])
    sort_command = 'perl %s --precise %s > %s' %(gff3sport_pl, gff_fn, out_fn)
    run_command(sort_command)
    return out_fn

def bed_file_split(bed_fn, num_chunks, tmp_path = TMP_PATH):
    '''Splits a bedfile into even chunks of files plus the reads. Each file having the same number of lines.
    Input:
    bed file name abspath
    chunks to split into
    tmp_path where the split files get saved into.
    Returns:
    List of abspath filenames for split bed files.'''
    infilename = bed_fn
    number_of_lines = 0
    with open(infilename, 'rb') as infile:
        for x in infile:
            number_of_lines = number_of_lines + 1
    print('in_file_line_number:', number_of_lines)
    chunk_size = number_of_lines // num_chunks
    print('target chunk_size:', chunk_size)
    files = []
    
    previous = 0
    for i,y in enumerate([x for x in range(0, number_of_lines, chunk_size)]):
        if y == 0:
            continue
                    
        else:
            temp_fn = os.path.join(TMP_PATH, '%s_%s' % (os.path.basename(bed_fn), i-1))
            with open(temp_fn, 'w') as temp_file:
                with open(infilename, 'r') as infile:
                    for line_number, line in enumerate(infile):
                        if line_number < y and line_number >= previous:
                            print(line.rstrip(), file=temp_file)
                        else:
                            continue
            previous = y
            files.append(temp_fn)
    if y < number_of_lines:
        temp_fn = os.path.join(TMP_PATH, '%s_%s' % (os.path.basename(bed_fn), i))
        with open(temp_fn, 'w') as temp_file:
                with open(infilename, 'r') as infile:
                    for line_number, line in enumerate(infile):
                        if line_number >= y:
                            print(line.rstrip(), file=temp_file)
                        else:
                            continue
    files.append(temp_fn)
    return files

def samcov_file_split(basefile_name, chunks, tmp_path = TMP_PATH):
    """This function takes the samcov file and generates the split filenames necessary to run the 
    parallelized run_sam_bedcov."""
    chunk_base_name = os.path.join(tmp_path, os.path.basename(basefile_name))
    return_list = []
    for x in range(0, chunks):
        return_list.append('%s_%s' % (chunk_base_name, x))
    return return_list

def combine_samcov_files(file_list):
    """This function concatenates samcov files from a tmp folder into the folder below. And removes all tmp files."""
    command = 'cat'
    combined_file = os.path.basename(file_list[0])[:-2]
    combined_file_fn = os.path.join(os.path.dirname(file_list[0]),'..', combined_file)
    for file in file_list:
        command = command + ' %s' % file
    command = command + ' > %s'  % combined_file_fn
    run_command(command)
    
def remove_files(file_list):
    for file in file_list:
        os.remove(file)

def run_sam_bedcov(bed_fn, bam_fn, samcov_fn, parallel = False):
    """This function runst samtools bedcov. It also has the option to do this in parallel spliting up the bed_cov file in equal junks plus one.
    Input:
        bed file name as Abspath.
        bam file name as Abspath.
        Samcov file as Abspath.
        Option to run thinks in parllel. Set to False. Provide the number of threads you want to run -1.
    Output:
        It generate the samcov files and saves them to file."""
    command = 'samtools bedcov %s %s > %s' % (bed_fn, bam_fn, samcov_fn)
    if parallel != False and type(parallel) == int:
        split_bed_files = bed_file_split(bed_fn, parallel)
        split_samcov_files = samcov_file_split(samcov_fn, len(split_bed_files), tmp_path=TMP_PATH)
        bam_files = [bam_fn]*len(split_bed_files)
        n_threads = len(split_bed_files)
        Parallel(n_jobs=n_threads)(delayed(run_sam_bedcov)(bed, bam, sam) for bed, bam, sam in zip(split_bed_files, bam_files, split_samcov_files))
        combine_samcov_files(split_samcov_files)
        remove_files(split_samcov_files)
        remove_files(split_bed_files)
    else:
        run_command(command)

def haploid_coverage(df):
    """This function tries to define the 'haploid coverage' of a samcov dataframe in case of two major coverage peaks.
    This should be run on a samcov file that is genereated by mapping againts primary and haplotigs. 
    In case of a single peak it will return its coverage. This functions needs some fixing most likely.
    It checks if the dataframe contains high coverage outliers and removes those for the analysis.
    Input:
        Samcov dataframe. 
    Output:
        Location of diploid coverage peak for diploids.
        Location of haploid coverage peak in case of haploids.
    """
    #function that first checks if the dataframe if trimmed down
    if df.ave_cov.max() > 10 * df.ave_cov.mean():
        print("Reducing dataframe to quantile 0.99")
        q_up = 0.99
        q_down = 1 - 0.99
        low = x.quantile(q_down) - 1.5*(x.quantile(q_up) - x.quantile(q_down))
        high = x.quantile(q_up) + 1.5*(x.quantile(q_up) - x.quantile(q_down))
        df = df[(df.ave_cov >= low) & (df.ave_cov <= high)]
    cov_bin = pd.cut(df.ave_cov, bins=200)
    value_counts = cov_bin.value_counts().sort_index()
    value_counts_df = value_counts.reset_index()
    local_max = argrelextrema(value_counts.values, np.greater, order=10, mode ='wrap')[0]
    sub_df = value_counts_df.loc[local_max]
    sub_df['mids'] = sub_df['index'].apply(lambda x: float(x.mid))
    sub_df = sub_df[sub_df.ave_cov > sub_df.ave_cov.max()*0.1]
    #now comes the heuristic
    print('These are the options for the peaks that are greater than 10% of the max peak.')
    print(sub_df)
    if len(sub_df.mids) == 1:
        #print("step 1")
        print("Please check if this is really the haploid coverage by looking at the 'haploid coverage' pcontigs when mapping against ph.")
        return 1, sub_df.mids[sub_df.mids.index[0]]
    elif len(sub_df.mids) == 2:
        if sub_df.mids[sub_df.mids.index[1]]/sub_df.mids[sub_df.mids.index[0]] > sub_df.mids[sub_df.mids.index[1]]:
        #print('step 2')
            return 2, sub_df.mids[sub_df.mids.index[1]]
        else:
            return 2, sub_df.mids[sub_df.mids.index[0]]
    elif len(sub_df.mids) > 2:
        #drop everything at the left edge
        mids = sub_df.mids[np.array(sub_df['mids']) >0.1*df.ave_cov.mean()]
        sub_df = sub_df[np.array(sub_df['mids']) >0.1*df.ave_cov.mean()]
        if len(mids) == 1:
            #print('something')
            return 1, mids[mids.index[0]]
            print("Please check if this is really the haploid coverage by looking at the 'haploid coverage' pcontigs when mapping against ph.")
            #continue
        elif len(mids) == 2:
            dvision = mids[mids.index[1]]/mids[mids.index[0]] 
            if 1.85 < dvision < 2.15:
                return 2, mids[mids.index[0]]
        elif len(mids) > 2:
            sub_df.reset_index(inplace=True, drop=True)
            max_index = sub_df[sub_df.ave_cov == sub_df.ave_cov.max()].index[0]
            #print(sub_df, max_index)
            return 3, mids[mids.index[max_index]]

def samcov_slurp(file_name, fil=True, low=0, high=1000000, quantile = False, norm = False, contig_fil = False):
    """Reads in a samcov file as a dataframe and can filter it by quantiles, high and low values.
    It can also normalize. Either by a given value or by the calculated haploid coverage if norm is set to true.
    If norm is set to false it will simply normalize by the overall mean coverage of the dataframe. This could be an issue
    for multi peak coverage plots.
    Input: 
        filename for samcov file.
    Options:
        fil: Filters the dataframe on low and high if quantile is False
        low: low filtering cut-off of coverage
        high: high filtering cut-off of coverage
            ONLY applicable if fil == True and quantile == False
        quantile: is the fraction used in the cut-off caclulcations.
        norm: False, True or integer. If False the norm coverage is simply calcuated dividing the ave cov by the mean of the ave_cov.
            If true normalization is done with the diploid coverage function. If int given norm is done with by dividing ave_cov with the integer.
        contig_fil: Can be either p or h to filter for pcontigs only or hcontigs only.
    """

    samcov_header = ['contig', 'start', 'stop', 'total_cov']
    df = pd.read_csv(file_name, sep='\t', header=None, names=samcov_header)
    df['ave_cov'] = df.total_cov/(df.stop-df.start)
    #rounder = pd.Series([0,0,0,0,2], index = df.columns)
    df.ave_cov = df.ave_cov.round()
    
    if fil == True and quantile == False:
        df = df[(df.ave_cov >= low) & (df.ave_cov <= high)]
    elif fil == True and type(quantile) == float:
        x = df.ave_cov
        q_up = quantile
        q_down = 1 - quantile
        low = x.quantile(q_down) - 1.5*(x.quantile(q_up) - x.quantile(q_down))
        high = x.quantile(q_up) + 1.5*(x.quantile(q_up) - x.quantile(q_down))
        df = df[(df.ave_cov >= low) & (df.ave_cov <= high)]
    if contig_fil == 'p':
        df = df[df.contig.str.startswith('pcontig')].copy()
    if contig_fil == 'h':
        df = df[df.contig.str.startswith('hcontig')].copy()
    if norm == False:
        df['norm_cov'] = df['ave_cov']/df['ave_cov'].mean()
    elif type(norm) == float or type(norm) == int:
        df['norm_cov'] = df['ave_cov']/float(norm)
    elif norm == True:
        peaks, hap_cov =  haploid_coverage(df)
        calcuated_coverage = float(np.round(hap_cov,3))
        df['norm_cov'] = df['ave_cov']/calcuated_coverage
        print("This is the calculated haploid coverage: %f" % calcuated_coverage)
        if peaks == 1:
            print("Please check if this is really the haploid coverage by looking at the 'haploid coverage' pcontigs when mapping against ph.")
    return df

def filter_on_norm_cov(samcov_df, cutoff=4):
    '''
    Filter a samcov_df on the norm cov column.
    Input: 
        samcov_df
        cufoff on which to filter the norm_cov column.
    Output:
        Filtered samcov_df.'''
    
    return samcov_df[samcov_df.norm_cov <=  cutoff].reset_index(drop=True)

def bin_value_filter_df(df, bins=500, fil=True):
    """Functions bins, a dataframe, valuecounts, finds mids, and filters bins with zero if fil == True.
    Inupt:
        Samcov dataframe.
        bins integer number for the number of bins.
        fil True/False for filtering out zero bins.
    Output:
        Dataframe with columns. Index, norm_cov and mids."""
    cov_bin = pd.cut(df.copy().norm_cov, bins=bins)
    value_counts = cov_bin.value_counts().sort_index()
    value_counts_df = value_counts.reset_index()
    value_counts_df['mids'] = value_counts_df['index'].apply(lambda x: float(x.mid))
    value_counts_df['norm_freq'] = value_counts_df['norm_cov']/ value_counts_df['norm_cov'].max()
    if fil == True:
        value_counts_df = value_counts_df[value_counts_df.norm_cov > 0]
    return value_counts_df.reset_index(drop=True)

def fill_plot_axis(xs, ys, ax, color, alpha=0.5):
    """Generates a plot axis object for a fillbetween plot."""
    ax.plot(xs, ys, marker='.', lw=2, color=color)
    d = scipy.zeros(len(xs))
    ax.fill_between(xs, ys, where=ys>=d, interpolate=True, color=color, alpha=alpha)
    return ax

def plot_mapping_annotation(ax, pallete, xstart_text, y_upper_lim):
    """Function that plots the legend text for each feature"""
    xincremend = xstart_text/20
    y_upper_spot = y_upper_lim/1.2
    yinremend = y_upper_lim/10
    ax[0,0].text(xstart_text, y_upper_lim/1.05, 'ref = primary contigs')
    ax[0,0].plot([xstart_text*1.1, xstart_text+(7*xincremend)], [y_upper_spot, y_upper_spot], lw =4, color = pallete[0])
    ax[0,0].text(xstart_text*0.97, y_upper_spot*0.98, 'p')


    #now plot the ph mapping with primeray mapping only
    ax[1,0].text(xstart_text, y_upper_lim/1.05, 'ref = primary contigs + haplotigs')
    ax[1,0].plot([xstart_text*1.1, xstart_text+(7*xincremend)], [y_upper_spot, y_upper_spot], lw =4, color = pallete[1])
    ax[1,0].text(xstart_text*0.97, y_upper_spot*0.98, 'p')
    ax[1,0].plot([xstart_text*1.1, xstart_text+(3*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = 'k')
    ax[1,0].plot([xstart_text+(5*xincremend), xstart_text+(7*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = 'k')
    ax[1,0].text(xstart_text*0.97, y_upper_spot*0.98*0.85, 'h')



    #now plot the ph mapping with haplotig mapping only


    ax[2,0].text(xstart_text, y_upper_lim/1.05, 'ref = primary contigs + haplotigs')
    ax[2,0].plot([xstart_text*1.1, xstart_text+(7*xincremend)], [y_upper_spot, y_upper_spot], lw =4, color = 'k')
    ax[2,0].text(xstart_text*0.97, y_upper_spot*0.98, 'p')
    ax[2,0].plot([xstart_text*1.1, xstart_text+(3*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = pallete[2])
    ax[2,0].plot([xstart_text+(5*xincremend), xstart_text+(7*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = pallete[2])
    ax[2,0].text(xstart_text*0.97, y_upper_spot*0.98*0.85, 'h')



    #now plot the ph mapping with primaries with overlaps

    ax[1,1].text(xstart_text, y_upper_lim/1.05, 'ref = primary contigs + haplotigs')
    ax[1,1].text(xstart_text, y_upper_lim/1.05, 'ref = primary contigs + haplotigs')
    ax[1,1].plot([xstart_text*1.1, xstart_text+(7*xincremend)], [y_upper_spot, y_upper_spot], lw =4, color = pallete[3])
    ax[1,1].plot([xstart_text+(3*xincremend), xstart_text+(5*xincremend)],[y_upper_spot, y_upper_spot] , lw =4, color = 'k')
    ax[1,1].text(xstart_text*0.97, y_upper_spot*0.98, 'p')
    ax[1,1].plot([xstart_text*1.1, xstart_text+(3*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = 'k')
    ax[1,1].plot([xstart_text+(5*xincremend), xstart_text+(7*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = 'k')
    ax[1,1].text(xstart_text*0.97, y_upper_spot*0.98*0.85, 'h')


    #now plot the ph mapping with primaries with overlaps


    ax[2,1].text(xstart_text, y_upper_lim/1.05, 'ref = primary contigs + haplotigs')
    ax[2,1].plot([xstart_text*1.1, xstart_text+(7*xincremend)], [y_upper_spot, y_upper_spot], lw =4, color = 'k')
    ax[2,1].plot([xstart_text+(3*xincremend), xstart_text+(5*xincremend)],[y_upper_spot, y_upper_spot] , lw =4, color = pallete[4])
    ax[2,1].text(xstart_text*0.97, y_upper_spot*0.98, 'p')
    ax[2,1].plot([xstart_text*1.1, xstart_text+(3*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = 'k')
    ax[2,1].plot([xstart_text+(5*xincremend), xstart_text+(7*xincremend)], [y_upper_spot*0.85, y_upper_spot*0.85], lw =4, color = 'k')
    ax[2,1].text(xstart_text*0.97, y_upper_spot*0.98*0.85, 'h')


    #add the labels
    fig.text(0.06, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    fig.text(0.5, 0.09, 'Normalized mapping Coverage',ha='center', va='center')

    #make plot [0,1] disappear
    ax[0,1].axis('off')
    text = '''
    Legend:\n
    ref == reference used for Illumina short read mapping
    p == primary contigs
    h == haplotigs
    color decodes genome regions that are analyzed for coverage
    e.g. pink indicates regions of primary contigs that have
    an overlapping haplotig
    '''
    ax[0,1].text(0, y_upper_lim/3, text)
    
    return ax
    
def df_to_bed_saved(df, lcutoff, samcov_fn):
    """Function that takes a dataframe, a lowcutoff for the normed coverage, and the samcov file.
    It filters the samcov on the lcutoff in the norm cov column. This is converted to a non-reduntand bed file using
    the bedtools merge function.
    Input: 
        samcov_df unfiltered but normalized samcov_df.
        lcutoff is a low cut off to filter the normed coverage on.
        samcov_fn for consistend naming should en on .samcov
    Ouput:
        fn_out of the final '.lowcovbed' file.
    The final output file should be filtered for low cov regions in the control sample."""
    fn = samcov_fn.replace('.samcov', '.lowcovbedtmp')
    df = df[df.norm_cov < lcutoff].reset_index()
    df.loc[:, ['contig', 'start', 'stop']].to_csv(fn, sep='\t', header=None, index=None)
    low_cov_bed = BedTool(fn)
    fn_out = samcov_fn.replace('.samcov', '.lowcov%sbed' % str(lcutoff).replace('.', ''))
    low_cov_bed.sort().merge().saveas(fn_out)
    os.remove(fn)
    return fn_out

def count_bases_in_bed(fn):
    """Function that sums the interval size of bedfile."""
    header = ['contig', 'start', 'stop']
    df = pd.read_csv(fn, sep='\t', header=None, names=header)
    df['interval'] = df.stop - df.start
    return df.interval.sum()