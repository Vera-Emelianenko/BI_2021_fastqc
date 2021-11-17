import collections
import csv
import math
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

input_file = "C:/Users/rusla/OneDrive/Desktop/1-5_S_F.fastq"
output_dir = "C:/Users/rusla/OneDrive/Desktop/seq/"


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
    plt.savefig(output_dir + 'Sequence_duplication_levels.png')
    fig_1

    ######
    # Write overrepresented sequences
    ######

    trunc_seq = []
    for seq in range(total_seq):
        t_seq = list_seq[seq][:50]
        trunc_seq.append(t_seq)
    dic_trunc_seq = collections.Counter(trunc_seq)

    with open(output_dir + 'Overrepresented_sequences.tsv', 'w',
              newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Sequence', 'Count', 'Percentage'])
        for key, values in dic_trunc_seq.items():
            if values * 100 / float(total_seq) >= 0.1:
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
    plt.savefig(output_dir + 'Adapter_content.png')
    fig_2


fastq_overseq(input_file)
