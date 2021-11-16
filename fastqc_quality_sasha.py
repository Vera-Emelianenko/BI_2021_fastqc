import matplotlib.pyplot as plt
#from numpy import histogram_bin_edges
import seaborn as sns
import pandas as pd

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

def per_sequence_quality_score():
    qualities = get_mean_quality()
    qual_dict = dict()
    qual_dict[1] = qualities
    df = pd.DataFrame(qual_dict)
    print(df)
    df.value_counts().sort_index().plot()
    plt.show()
#per_sequence_quality_score()

def get_quality_base():
    input = get_seq_quality()
    bases = [[] for j in range(0, len(input[1]))]
    one_col_length = len(input)*len(input[1])
    bases_one_col = [i for i in range(1, one_col_length+1)]
    for i in range(0, len(input)):
        for j in range(0, len(input[1])):
            bases[j].append(input[i][j])
            bases_one_col.append(input[i][j])
    return(bases)
print(get_quality_base())

def get_quality_base_one_col():
    input = get_seq_quality()
    bases_one_col = {"quality": [], "read_number": [], "position": []}
    for i in range(0, len(input)):
        for j in range(0, len(input[i])):
            bases_one_col["quality"].append(input[i][j])
            bases_one_col["read_number"].append(i)
            bases_one_col["position"].append(j)
    return(bases_one_col)
print(get_quality_base_one_col())

def per_base_sequence_quality():
    df = pd.DataFrame(get_quality_base_one_col())
    print(df)
    sns.boxplot(x = "position", y = "quality", data = df)
    plt.show()
per_base_sequence_quality()

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
#per_tile_quality()