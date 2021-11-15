import matplotlib.pyplot as plt
#from numpy import histogram_bin_edges
import seaborn as sns
import pandas as pd

#def get_seq_quality(input):
def get_seq_quality():
    #with open(input) as input_file:
        #input_fastq_list = input_file.readlines()
    input_fastq_list = ['@SEQ_ID', 'ATGCTGC', '+', '@AB!CDA@', '@SEQ_ID2', 'ATGCGGC', '+', '@ABDC>A@', '@SEQ_ID3', 'AAACTGC', '+', '@A>>CDA@', 
    '@SEQ_ID4', 'AAACTGC', '+', '@A>>CDA@']
    seq_quality = []
    for i in range(0, len(input_fastq_list)):
        if (i+1)%4 == 0:
            decoded = []
            for el in input_fastq_list[i]:
                decoded.append(ord(el))
            seq_quality.append(decoded)
    return(seq_quality)
print(get_seq_quality())

#don't forget to include input!!!
def get_mean_quality(input = get_seq_quality()):
    mean_quality = []
    for i in input:
        mean_quality.append(sum(i)/len(i))
    return(mean_quality)
print(get_mean_quality())

#def per_sequence_quality_score():
    #qualities = get_mean_quality()
    #sns.histplot(qualities, bins=1, binrange=(10,99), fill=False, color ='red')
    #plt.show()
#per_sequence_quality_score()

def get_quality_base():
    input = get_seq_quality()
    bases = [[] for j in range(0, len(input[1]))]
    for i in range(0, len(input)):
        for j in range(0, len(input[1])):
            bases[j].append(input[i][j])
    return(bases)
print(get_quality_base())

def get_tile():
    input_fastq_list = ['@SIM:1:FCX:1:1:6329:1045', 'ATGCTGC', '+', '@AB!CDA@', '@@SIM:1:FCX:1:1:6329:1045:GATTACT', 'ATGCGGC', '+', 
    '@ABDC>A@', '@SIM:1:FCX:1:3:6329:1045:GATTACT', 'AAACTGC', '+', '@A>>CDA@', 
    '@SIM:1:FCX:1:3:6329:1045:GATTACT', 'AAACTGC', '+', '@A>>CDA@']
    tiles = []
    for i in range(0, len(input_fastq_list)):
        if (i)%4 == 0:
            tiles.append(input_fastq_list[i].split(':')[4])
    return(tiles)
print(get_tile())

def base_quality_per_tile():
    tiles = get_tile()
    tiles_dict_quality = dict.fromkeys(tiles)
    tiles_count = dict.fromkeys(tiles)
    seq = get_seq_quality()
    for i in tiles:
        tiles_count[i] = tiles.count(i)
        if tiles_dict_quality[i] is None:
           tiles_dict_quality[i] =  seq[tiles.index(i)]
        else:
            tiles_dict_quality[i] = [sum(x) for x in zip(tiles_dict_quality[i], seq[tiles.index(i)])]
        print(tiles_dict_quality)
    for key in tiles_dict_quality:
        for item in tiles_dict_quality[key]: 
            ind = tiles_dict_quality[key].index(item)
            item = item/(tiles_count[key])
            tiles_dict_quality[key][ind] = item
    return(tiles_dict_quality)
print(base_quality_per_tile())





def per_tile_quality():
    df = pd.DataFrame.from_dict(base_quality_per_tile(), orient='index')
    df = df.set_axis([i for i in range(1, len(df.columns)+1)], axis='columns')
    df = df.reindex(index=df.index[::-1])
    print(df)
    sns.heatmap(df)
    sns.color_palette("coolwarm", as_cmap=True)
    plt.show()
per_tile_quality()