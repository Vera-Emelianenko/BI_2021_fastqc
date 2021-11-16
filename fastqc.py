import os
import argparse
from Bio import SeqIO
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
# for testing purposes only
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
        temp_list.append([fastq_read.id, str(fastq_read.seq)])
    df = pd.DataFrame(temp_list, columns=["id", "seq"])
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

    fig.savefig(os.path.join(args.outdir, "per_sequence_gc_content.png"), format='png', dpi=600)


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

    fig.savefig(os.path.join(args.outdir, "per_base_sequence_content.png"), format='png', dpi=600)


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

    fig.savefig(os.path.join(args.outdir, "per_base_n_content.png"), format='png', dpi=600)


def main():
    start_time = time()
    df = fastq_to_dataframe(args.input)
    plot_gc_content(df)
    plot_sequence_content(df)
    plot_n_content(df)
    end_time = time()
    seconds_elapsed = end_time - start_time
    hours, rest = divmod(seconds_elapsed, 3600)
    minutes, seconds = divmod(rest, 60)
    print(f'''Analysis complete for {args.input} in {hours} hours {minutes}
          minutes {seconds} seconds. Results written to "{args.outdir}/"''')


main()
