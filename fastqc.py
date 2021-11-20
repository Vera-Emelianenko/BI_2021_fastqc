import os
import argparse
from Bio import SeqIO
import pandas as pd
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.patches as patches
import matplotlib.ticker as plticker
from matplotlib.ticker import (AutoLocator, AutoMinorLocator, MultipleLocator)
import seaborn as sns
from time import time



parser = argparse.ArgumentParser(description="FastQC analog as a homework project")
parser.add_argument("-o", "--outdir",
                    help="Specify directory in which output has to be created,  default ./", default=".")
parser.add_argument("-i", "--input", help="Reads data in fastq format")
args, unknown = parser.parse_known_args()


def fastq_to_dataframe(filename):

    # check if extension is right, then parse the file

    ext = os.path.splitext(filename)[1]
    if ext == '.fastq':
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    else:
        return("The programm takes only fastq files! Try again.")

    # create a temporary list with lists of two elements: id and sequence
    # then create a dataframe and add the column "length"

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

    ax.hist(df.gc_content, bins=40, color="red", histtype="step", lw=2, label="GC count per read")
    ax.set_xlim((0, 100))
    ax.set_xlabel("Mean GC content (%)", fontsize="16")
    ax.tick_params(axis='y', labelcolor="black")
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    x = np.arange(1, 100, 0.1)
    norm_fig = [norm_gc(i, df.gc_content.mean(), df.gc_content.std()) for i in x]
    ax2 = ax.twinx()  # the second y axis
    ax2.plot(x, norm_fig, label="Theoretical distribution", color="blue", lw=3)
    ax2.set_ylim(bottom=0)
    ax2.axes.get_yaxis().set_visible(False)

    ax.set_title("GC distribution over all sequences", size=20)
    ax.legend(loc=2, bbox_to_anchor=(0.69, 1.0))
    ax2.legend(loc=2, bbox_to_anchor=(0.69, 0.95))

    fig.savefig(os.path.join(args.outdir, os.path.splitext(args.input)[0] + "_per_sequence_gc_content.png"),
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

    fig.savefig(os.path.join(args.outdir, os.path.splitext(args.input)[0] + "_per_base_sequence_content.png"),
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

    fig.savefig(os.path.join(args.outdir, os.path.splitext(args.input)[0] + "_per_base_n_content.png"),
                format='png', dpi=300)



def plot_length_distribution(df):

    fig, ax = plt.subplots(figsize=(15, 10))

    plt.rcParams['font.size'] = '16'
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    m1 = min(df.length)
    m2 = max(df.length)

    ax.hist(df.length, bins=np.arange(m1, m2+1)-0.5, rwidth=0.5, color="red", lw=2, label="Sequence length")
    ax.set_xlim((m1, m2+1))
    ax.set_xlabel("Sequence length (bp)", fontsize="16")
    ax.tick_params(axis='y', labelcolor="black")
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    ax.set_title("Distribution of sequence lengths over all sequences", size=20)
    ax.legend(loc=2)

    fig.savefig(os.path.join(args.outdir, os.path.splitext(args.input)[0] + "_seq_length_distribution.png"),
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
    import math
    mpl.rcParams.update(mpl.rcParamsDefault)
    qualities = get_mean_quality(input_fastq_list)
    df = pd.DataFrame({"quality": qualities}, columns=["quality"])
    quality_counts = pd.DataFrame(df["quality"].value_counts().sort_index(), columns=["quality"])
    fg, ax = plt.subplots()
    ax.plot(quality_counts, label="Average quality per read")
    ax.xaxis.set_major_locator(AutoLocator())
    ax.yaxis.set_major_locator(AutoLocator())
    plt.title('Quality score distribution over all sequences')
    plt.xlabel("Quality")
    plt.legend(loc='upper right')
    
    plt.savefig(os.path.join(args.outdir, "per_sequence_quality_score.png"), format='png', dpi=1000)


# Per base quality is extracted, put to list with lists for each read position (bases)
def get_quality_base(input_fastq_list):
    input = get_seq_quality(input_fastq_list)
    bases = [[] for j in range(0, len(input[1]))]
    for i in range(0, len(input)):
        for j in range(0, len(input[1])):
            bases[j].append(input[i][j])
    return(bases)


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
    df = pd.DataFrame(get_quality_base_one_col(), columns=["quality", "position"])
    sns.set_style("whitegrid")
    qual_mean = pd.DataFrame(df.groupby('position').mean().reset_index(), columns=["quality", "position"])
    per_base_plot = sns.boxplot(x="position", y="quality", data=df, hue = "position", showfliers = False, width=80)
    # Areas of background color are defined depending on y axis values
    per_base_plot.axhspan(0, 20, color='#EF9A9A', alpha=0.4)
    per_base_plot.axhspan(20, 28, color=(0.9, 1, 0.5, 1), alpha=0.4)
    per_base_plot.axhspan(28, 34, color='#388E3C', alpha=0.4)

    per_base_plot.set(title='Quality scores across all bases (Illumina>v1.3 encoding)', xlabel='Position in read (bp)')
    per_base_plot.xaxis.set_major_locator(AutoLocator())
    per_base_plot.get_legend().remove()
    plt.plot(qual_mean["quality"])
    plt.ylim(0, 40)
    plt.xlim(0, 40)
    
    plt.savefig(os.path.join(args.outdir, "per_base_sequence_quality.png"), format='png', dpi=1000)


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
    df = pd.DataFrame.from_dict(base_quality_per_tile(input_fastq_list), orient='index')
    df = df.set_axis([i for i in range(1, len(df.columns)+1)], axis='columns')
    df = df.reindex(index=df.index[::-1])
    per_tile_plot = sns.heatmap(df, cmap="RdBu")
    per_tile_plot.xaxis.set_major_locator(AutoLocator())
    per_tile_plot.set(title='Quality per tile', xlabel='Position in read (bp)', ylabel='Tile')
    
    plt.savefig(os.path.join(args.outdir, "per_tile_quality.png"), format='png', dpi=1000)


def main():
    start_time = time()
    df = fastq_to_dataframe(args.input)
    plot_gc_content(df)
    plot_sequence_content(df)
    plot_n_content(df)
    plot_length_distribution(df)
    end_time = time()
    seconds_elapsed = end_time - start_time
    hours, rest = divmod(seconds_elapsed, 3600)
    minutes, seconds = divmod(rest, 60)
    print(f'''Analysis complete for {args.input} in {hours} hours {minutes}
          minutes {seconds} seconds. Results written to "{args.outdir}/"''')

    with open(args.input) as input_fastq:
        fastq_list = input_fastq.read().splitlines()
    get_seq_quality(fastq_list)
    per_sequence_quality_score(fastq_list)
    per_base_sequence_quality(fastq_list)
    try:
        get_tile(fastq_list)
        per_tile_quality(fastq_list)
    except IndexError:
        print('No tile info provided in fastq file, unable to generate per tile plot')
    fastq_overseq(args.input)


main()
