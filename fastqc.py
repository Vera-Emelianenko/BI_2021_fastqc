import os
import sys
import argparse
from Bio import SeqIO
import pandas as pd
import math
import collections
import csv
import numpy as np
from scipy.signal import find_peaks as find_peaks
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoLocator, AutoMinorLocator, MaxNLocator)
import seaborn as sns
from time import time

if len(sys.argv) == 1:
    sys.exit("No arguments provided. Print python fastqc.py -h to see help message")
parser = argparse.ArgumentParser(description="FastQC analog as a homework project")
parser.add_argument("-o", "--outdir",
                    help="Specify directory in which output has to be created,  default ./", default=".")
parser.add_argument("-i", "--input", help="Reads data in fastq format")
args, unknown = parser.parse_known_args()


def check_path(input_path, output_path):
    if not os.path.exists(input_path):
        sys.exit("Input file not found")
    elif not os.path.exists(output_path):
        sys.exit("Output directory not found")
    else:
        print(f'Started analysis of {args.input}')


def fastq_to_dataframe(filename):

    # check if extension is right, then parse the file
    ext = os.path.splitext(filename)[1]
    if ext == '.fastq':
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    else:
        sys.exit("The programm takes only .fastq files! Try again")

    temp_list = []
    for fastq_read in fastq_parser:
        temp_list.append(str(fastq_read.seq))
    df = pd.DataFrame(temp_list, columns=["seq"])
    df["length"] = df.seq.str.len()

    def gc_content(seq):
        gc_content = round(((seq.count('C') + seq.count('G'))/len(seq)*100), 2)
        return gc_content
    df["gc_content"] = df.seq.apply(gc_content)

    return df


def norm_gc(x, mean, sd):
    """create a normal distribution at given GC-count of the organism"""

    variance = float(sd)**2
    denominator = (2*math.pi*variance)**0.5
    numerator = math.exp(-(float(x)-float(mean))**2/(2*variance))
    return numerator/denominator


def plot_gc_content(df):

    fig, ax = plt.subplots(figsize=(15, 10))

    plt.rcParams['font.size'] = '16'
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    gc = np.array(df.gc_content)-0.025
    hist, bin_edges = np.histogram(gc, 40)
    bin_edges = bin_edges[1:]
    peaks, _ = find_peaks(hist)

    ax.plot(bin_edges, hist, color="red", lw=2, label="GC count per read")
    ax.set_xlim((0, 100))
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Mean GC content (%)", fontsize="16")
    ax.tick_params(axis='y', labelcolor="black")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    x = np.arange(1, 100, 0.1)
    norm_fig = [norm_gc(i, df.gc_content.mean(), df.gc_content.std()) for i in x]
    ax2 = ax.twinx()  # the second y axis
    ax2.plot(x, norm_fig, label="Theoretical distribution", color="blue", lw=2)
    ax2.set_ylim(bottom=0)
    ax2.axes.get_yaxis().set_visible(False)

    ax.set_title("GC distribution over all sequences", size=20)
    ax.legend(loc=2, bbox_to_anchor=(0.69, 1.0))
    ax2.legend(loc=2, bbox_to_anchor=(0.69, 0.95))

    fig.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_per_sequence_gc_content.png"),
                format='png', dpi=300)


def base_content(df):

    base_dict = dict()

    for seq in df.seq:
        j = 1
        for base in seq:
            if j not in base_dict.keys():
                base_dict[j] = list()
            base_dict[j].append(base)
            j += 1

    T_count = list()
    C_count = list()
    A_count = list()
    G_count = list()
    N_count = list()

    for value in base_dict.values():
        freq_t = 100*value.count("T")/len(value)
        freq_c = 100*value.count("C")/len(value)
        freq_a = 100*value.count("A")/len(value)
        freq_g = 100*value.count("G")/len(value)
        freq_n = 100*value.count("N")/len(value)
        T_count.append(freq_t)
        C_count.append(freq_c)
        A_count.append(freq_a)
        G_count.append(freq_g)
        N_count.append(freq_n)

    count_list = [T_count, C_count, A_count, G_count, N_count]

    return count_list


def plot_sequence_content(df):

    count_list = base_content(df)

    fig, ax = plt.subplots(figsize=(15, 10))

    plt.rcParams['font.size'] = '16'
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    ax.plot(range(1, len(count_list[0])+1), count_list[0], color='red', label="%T")
    ax.plot(range(1, len(count_list[1])+1), count_list[1], color='blue', label="%C")
    ax.plot(range(1, len(count_list[2])+1), count_list[2], color='green', label="%A")
    ax.plot(range(1, len(count_list[3])+1), count_list[3], color='black', label="%G")

    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    ax.set_ylim((0, 100))
    ax.legend(loc='upper right', labelcolor='linecolor')

    ax.set_xlabel('Position in read (bp)')
    ax.set_title("Sequence content across all bases", size=20)

    fig.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_per_base_sequence_content.png"),
                format='png', dpi=300)


def plot_n_content(df):

    count_list = base_content(df)

    fig, ax = plt.subplots(figsize=(15, 10))

    plt.rcParams['font.size'] = '16'
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    ax.plot(range(1, len(count_list[4])+1), count_list[4], 'red', label="%N", lw=3)

    ax.set_ylim((0, 100))
    ax.legend(loc='upper right')

    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    ax.set_xlabel('Position in read (bp)')
    ax.set_title("N content across all bases", size=20)

    fig.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_per_base_n_content.png"),
                format='png', dpi=300)


def plot_length_distribution(df):

    fig, ax = plt.subplots(figsize=(15, 10))

    plt.rcParams['font.size'] = '16'
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    len_distr = np.array(df.length)-0.025
    hist, bin_edges = np.histogram(len_distr, 40)
    bin_edges = bin_edges[1:]
    peaks, _ = find_peaks(hist)

    ax.plot(bin_edges, hist, color="red", lw=2, label="Sequence length")
    ax.set_xlim((min(df.length)-1, max(df.length)+1))
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Sequence length (bp)", fontsize="16")
    ax.tick_params(axis='y', labelcolor="black")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    ax.set_title("Distribution of sequence lengths over all sequences", size=20)
    ax.legend(loc=2)

    fig.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_seq_length_distribution.png"),
                format='png', dpi=300)


# creating small internal test

internal_test = ['@SIM:1:FCX:1:1:6329:1045:GATTACT', 'ATGCTGC', '+', '@AB!CDA@']
internal_test += ['@@SIM:1:FCX:1:1:6329:1045:GATTACT', 'ATGCGGC', '+', '@ABDC>A@']
internal_test += ['@SIM:1:FCX:1:3:6329:1045:GATTACT', 'AAACTGC', '+', '@A>>CDA@']
internal_test += ['@SIM:1:FCX:1:3:6329:1045:GATTACT', 'AAACTGC', '+', '@A>>CDA@']


# The fastq lines corresponding to quality are extracted, decoded to scores and
# put to list with lists each of which corresponds to 1 read
def get_seq_quality(input_fastq_list):
    seq_quality = []
    for i in range(0, len(input_fastq_list)):
        if (i+1) % 4 == 0:
            decoded = []
            for el in input_fastq_list[i]:
                decoded.append(ord(el)-33)
            seq_quality.append(decoded)
    return(seq_quality)


# Mean quality per read is calculated and put to list
def get_mean_quality(input_fastq_list):
    input = get_seq_quality(input_fastq_list)
    mean_quality = []
    for i in input:
        mean_quality.append(sum(i) / len(i))
    return(mean_quality)


# Number of quality per read (mean) values occurences is calculated, put to data frame and plotted
def per_sequence_quality_score(input_fastq_list):
    mpl.rcParams.update(mpl.rcParamsDefault)
    qualities = get_mean_quality(input_fastq_list)
    qual_array = np.array(qualities)
    hist, bin_edges = np.histogram(qual_array, 10)
    bin_edges = bin_edges[1:]
    plt.plot(bin_edges, hist, lw=2, label="Average quality per read")
    plt.title('Quality score distribution over all sequences')
    plt.xlabel("Quality")
    plt.legend(loc='upper right')
    plt.show()

    plt.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_per_sequence_quality_score.png"),
                format='png', dpi=300)


# Per base quality is extracted and put to long-format column dictionary (bases_one_col)
def get_quality_base_one_col(input_fastq_list):
    input = get_seq_quality(input_fastq_list)
    bases_one_col = {"quality": [], "read_number": [], "position": []}
    for i in range(0, len(input)):
        for j in range(0, len(input[i])):
            bases_one_col["quality"].append(input[i][j])
            bases_one_col["read_number"].append(i)
            bases_one_col["position"].append(j + 1)
    return(bases_one_col)


# Per base quality is plotted, line representing of per base qualities (qual_mean) is plotted
def per_base_sequence_quality(input_fastq_list):
    mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set_style("whitegrid")
    df = np.array(list(get_quality_base_one_col(input_fastq_list).values()))
    df = np.transpose(df)
    df = df[df[:, 0].argsort()]
    df2 = list(np.split(df[:, 1], np.unique(df[:, 0], return_index=True)[1][1:]))
    qual_mean = []
    for x in range(len(df2)):
        qual_mean.append(np.mean(df2[x]))
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(15, 10))
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    per_base_plot = sns.boxplot(x=df[:, 0],  y=df[:, 1], showfliers=False,
                                width=0.5, color='yellow', ax=ax)
    # Areas of background color are defined depending on y axis values
    per_base_plot.axhspan(0, 20, color='#EF9A9A', alpha=0.4)
    per_base_plot.axhspan(20, 28, color=(0.9, 1, 0.5, 1), alpha=0.4)
    per_base_plot.axhspan(28, 34, color='#388E3C', alpha=0.4)
    per_base_plot.axhspan(28, 40, color='#388E3C', alpha=0.4)

    per_base_plot.set(title='Quality scores across all bases (Illumina>v1.3 encoding)', xlabel='Position in read (bp)')
    ax.set_title('Quality scores across all bases (Illumina>v1.3 encoding)', size=20)
    ax.set_xlabel('Position in read (bp)', fontsize=16)
    ax.set_ylabel('Quality', fontsize=16)
    per_base_plot.xaxis.set_major_locator(AutoLocator())
    per_base_plot.get_legend().remove()

    plt.plot(qual_mean)
    plt.ylim(0, 40)
    plt.xlim(0, 40)

    plt.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_per_base_sequence_quality.png"),
                format='png', dpi=300)


# The tile id is extracted as 5th column in 1st line for each fastq entry
def get_tile(input_fastq_list):
    tiles = []
    for i in range(0, len(input_fastq_list)):
        if (i) % 4 == 0:
            tiles.append(input_fastq_list[i].split(':')[4])
    return(tiles)


# Per tile quality is extracted, quality values [items] are put to dictionary according to tiles ids [keys]
def base_quality_per_tile(input_fastq_list):
    tiles = get_tile(input_fastq_list)
    tiles_dict_quality = dict.fromkeys(tiles)
    tiles_count = dict.fromkeys(tiles)
    seq = get_seq_quality(input_fastq_list)
    for i in tiles:
        tiles_count[i] = tiles.count(i)
        if tiles_dict_quality[i] is None:
            tiles_dict_quality[i] = seq[tiles.index(i)]
        else:
            tiles_dict_quality[i] = [sum(x) for x in zip(tiles_dict_quality[i], seq[tiles.index(i)])]
    for key in tiles_dict_quality:
        for item in tiles_dict_quality[key]:
            ind = tiles_dict_quality[key].index(item)
            item = item/(tiles_count[key])
            tiles_dict_quality[key][ind] = item
    return(tiles_dict_quality)


# Per tile quality is plotted as heatmap
def per_tile_quality(input_fastq_list):
    mpl.rcParams.update(mpl.rcParamsDefault)
    qual_array = np.array(list(base_quality_per_tile(input_fastq_list).values()))
    plt.imshow(qual_array, cmap='hot', interpolation='nearest')
    plt.title('Quality per tile')
    plt.xlabel("Position in read (bp)")
    plt.ylabel("Tile")
    plt.show()

    plt.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_per_tile_quality.png"),
                format='png', dpi=300)


def fastq_overseq(input_file):

    ######
    # Read input file
    ######

    list_seq = []
    for record in SeqIO.parse(input_file, "fastq"):
        list_seq.append(record.seq)
    dic_seq = collections.Counter(list_seq)
    total_seq = len(list_seq)
    uniq_seq = len(dic_seq)
    freq_list = list(dic_seq.values())

    ######
    # Plot sequence duplication levels
    ######

    def between(list1, low, high):
        list2 = []
        for i in list1:
            if(i >= low and i < high):
                list2.append(i)
        return list2

    counting_borders = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000,
                        5000, 10000, math.inf]
    y_blue = []
    y_red = []
    x_values = []
    for i in range(len(counting_borders) - 1):
        low = counting_borders[i]
        high = counting_borders[i+1]
        freq_freq = between(freq_list, low, high)
        freq_blue = sum(freq_freq) / float(total_seq)
        freq_red = len(freq_freq) / float(uniq_seq)
        y_blue.append(freq_blue)
        y_red.append(freq_red)
        x_values.append(i)

    labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '>10', '>50',
              '>100', '>500', '>1k', '>5k', '>10k']
    fig_1 = plt.figure(figsize=(10, 7))
    plt.plot(x_values, np.array(y_blue)*100, label='% Total sequences')
    plt.plot(x_values, np.array(y_red)*100, label='% Deduplicated sequences')
    plt.xticks(ticks=x_values, labels=labels)
    plt.ylim(0, 100)
    plt.grid(linewidth=0.5)
    plt.title('Percent of seqs remaining if deduplicated ' +
              str(round(uniq_seq / float(total_seq) * 100, 2)) + ' %')
    plt.legend()
    plt.show()
    plt.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + "_sequence_duplication_levels.png"),
                format='png', dpi=300)
    fig_1

    ######
    # Write overrepresented sequences
    ######

    trunc_seq = []
    for seq in range(total_seq):
        t_seq = list_seq[seq][:50]
        trunc_seq.append(t_seq)
    dic_trunc_seq = collections.Counter(trunc_seq)

    with open(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + '_overrepresented_sequences.tsv'), 'w',
              newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Sequence', 'Count', 'Percentage'])
        for key, values in dic_trunc_seq.items():
            if values * 100 / float(total_seq) >= 0.1 and values > 1:
                tsv_writer.writerow([key, values,
                                     values * 100 / float(total_seq)])

    ######
    # Plot adapter content
    ######

    ad_dict = {'llumina Universal Adapter': 'AGATCGGAAGAG',
               'llumina Small RNA Adapter': 'ATGGAATTCTCG',
               'Nextera Transposase Sequence': 'CTGTCTCTTATA',
               'SOLID Small RNA Adapter': 'CGCCTTGGCCGT'}
    lens_seq = []
    for i in range(len(list_seq)):
        len_seq = len(list_seq[i])
        lens_seq.append(len_seq)
    max_len_seq = max(lens_seq)
    fig_2 = plt.figure(figsize=(10, 7))
    for ad_name, ad_seq in ad_dict.items():
        pos_seq = [0] * (max_len_seq - len(ad_seq))
        for seq in list_seq:
            for nucl in range(len(seq) - len(ad_seq)):
                substr = seq[nucl:nucl + len(ad_seq)]
                if substr == ad_seq:
                    pos_seq[nucl] += 1
        pos_seq_cum = [pos_seq[0]]
        for i in range(len(pos_seq) - 1):
            pos_seq_cum.append(pos_seq_cum[-1] + pos_seq[i + 1])
        perc_pos_seq_cum = np.array(pos_seq_cum) * 100 / float(total_seq)
        plt.plot(range(1, (max_len_seq - len(ad_seq) + 1)), perc_pos_seq_cum,
                 label=ad_name)
    plt.ylim(0, 100)
    plt.grid(linewidth=0.5)
    plt.title('% Adapter')
    plt.legend()
    plt.show()
    plt.savefig(os.path.join(args.outdir, os.path.basename(args.input)[:-6] + '_adapter_content.png'),
                format='png', dpi=300)
    fig_2


def print_end_time(start_time):
    end_time = time()
    seconds_elapsed = end_time - start_time
    hours, rest = divmod(seconds_elapsed, 3600)
    minutes, seconds = divmod(rest, 60)
    print(f'Analysis completed for {args.input}. Results written to {args.outdir}')
    print(f'Time: {hours} hours {minutes} minutes {seconds} seconds. ''')


def print_base_statistics(df, input_file, output_dir):

    '''writes down a .tsv file with basic statistics'''

    filename = os.path.basename(input_file)
    average_gc_content = str(round(df.gc_content.mean(), 1))
    min_length = str(df.length.min())
    max_length = str(df.length.max())
    mean_length = str(round(df.length.mean(), 1))
    if max_length == min_length:
        sequence_length_range = max_length
    else:
        sequence_length_range = min_length + '-' + max_length
    number_of_sequences = str(len(df))
    basic_statistics_dict = {
        'Filename': filename,
        'Total Sequences': number_of_sequences,
        'Mean sequence length': mean_length,
        'Sequence length range': sequence_length_range,
        'Average GC-content, %': average_gc_content
    }

    with open(os.path.join(output_dir, os.path.basename(filename)[:-6] + '_basic_statistics.tsv'), 'w',
              newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Measure', 'Value'])
        for keys, values in basic_statistics_dict.items():
            tsv_writer.writerow([keys, values])
    return None


def main():
    start_time = time()
    check_path(args.input, args.outdir)
    df = fastq_to_dataframe(args.input)
    print('Plotting GC content...')
    plot_gc_content(df)
    print('Plotting sequence content...')
    plot_sequence_content(df)
    plot_n_content(df)
    print('Plotting length distribution...')
    plot_length_distribution(df)
    print('Plotting sequence quality...')
    with open(args.input) as input_fastq:
        fastq_list = input_fastq.read().splitlines()
    get_seq_quality(fastq_list)
    per_sequence_quality_score(fastq_list)
    per_base_sequence_quality(fastq_list)
    try:
        get_tile(fastq_list)
        plt.close("all")
        per_tile_quality(fastq_list)
    except IndexError:
        print('No tile info provided in fastq file, unable to generate per tile plot')
    print('Calculating overrepresented sequences...')
    fastq_overseq(args.input)
    print_base_statistics(df, args.input, args.outdir)
    print_end_time(start_time)


main()
