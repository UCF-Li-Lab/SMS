import pysam
import os.path
import argparse
import Constants
import sys
import numpy as np
def read_polymorphicsSite(inName):
    data = []
    with open('%s' % inName) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            line = [int(one) for one in line]
            data.append(line)
    return data

import pandas as pd

# def read_polymorphicsSite(filename):
#     return pd.read_csv(filename,index_col=0,names=['A','C','G','T'],sep='\t')

# def read_polymorphicsSite(inName):
#
#     df = reader(inName)
#     # for index, row in df.iterrows():
#
#
#
#     return


# filter position with total reads less than 10% of all position.
# def filter(data):
#     res = []
#     # filter the position with total number of reads less than 10%
#     total_read = [sum(item[1:]) for item in data]
#     quantile = np.quantile(total_read, 0.1)
#     total_read.sort()
#
#     for item in data:
#         if sum(item[1:]) <= quantile:
#             continue
#         else:
#             res.append(item[:])
#
#     return res

def filter(data):
    res = []
    # filter the position with total number of reads less than 10%
    total_read = [sum(item[1:]) for item in data]
    ten_percentage = 0.1*sum(total_read)/len(total_read)

    for item in data:
        if sum(item[1:]) <= ten_percentage:
            continue
        else:
            res.append(item[:])
    return res

def updated_genome(polyCounts, genomeFileLoc):
    genome_name = []
    genome = ''
    with open('%s' % genomeFileLoc) as f:
        for line in f:
            if line.startswith('>'):
                genome_name.append(line.strip())
                pass
            else:
                genome += line.strip()

    genome_dict = {}
    for i in range(len(genome)):
        genome_dict[i] = genome[i]

    nucl_order = ['A','C','G','T']
    for item in polyCounts:
        pos = item[0]
        A_fre = item[1]
        C_fre = item[2]
        G_fre = item[3]
        T_fre = item[4]
        temp_sum = sum([int(i) for i in item[1:]])
        if temp_sum == int(A_fre):
            genome_dict[pos]='A'
        if temp_sum == int(C_fre):
            genome_dict[pos]='C'
        if temp_sum == int(G_fre):
            genome_dict[pos]='G'
        if temp_sum == int(T_fre):
            genome_dict[pos]='T'

        max_value = max([int(i) for i in item[1:]])

        if temp_sum > max_value:
            ref_nul = genome_dict[pos]
            index1 = nucl_order.index(ref_nul)
            temp_fre = int(item[1:][index1])

            if temp_fre == 0:
                index2 = item[1:].index(max_value)
                nul_temp = nucl_order[index2]
                genome_dict[pos] = nul_temp

    f1 = open(genomeFileLoc + '_updated','w')
    val_str = ''
    for key,val in genome_dict.items():
        val_str +=val
    f1.write(genome_name[0] + '\n')
    for i in range(0, len(val_str), 80):
        temp_seq = val_str[i:i + 80]
        f1.write(temp_seq + '\n')
    f1.close()

def filter_four_zero(res):
    polymorphicSites = []
    for loc in res:
        flag_nonZero = 0
        for nucl in loc[1:]:
            nucl_count = nucl
            if nucl_count != 0:
                flag_nonZero += 1
        if flag_nonZero > 1:
            polymorphicSites.append(loc)
    return polymorphicSites

def getMergedPolymorphicSites(genome_len,reads_file_list, res_dir,genomeFileLoc):
    res = []
    for _ in range(genome_len):
        res.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})

    res_count = []
    for _ in range(genome_len):
        res_count.append({'A': [], 'C': [], 'G': [], 'T': []})

    # res_matrix = {}
    # str1 = []
    for reads_file in reads_file_list:
        # generate polymorphic sites
        print(reads_file)
        sub_polymorphicSites = read_polymorphicsSite(reads_file)
        for item in sub_polymorphicSites:
            pos = item[0]
            A_count = item[1]
            C_count = item[2]
            G_count = item[3]
            T_count = item[4]

            if A_count > 0:
                res[pos]['A'] += A_count
                res_count[pos]['A'].append(A_count)
            if C_count > 0:
                res[pos]['C'] += C_count
                res_count[pos]['C'].append(C_count)
            if G_count > 0:
                res[pos]['G'] += G_count
                res_count[pos]['G'].append(G_count)
            if T_count > 0:
                res[pos]['T'] += T_count
                res_count[pos]['T'].append(T_count)

            # temp_str = [pos,[reads_file_list.index(reads_file),A_count,C_count,G_count,T_count]]
            # str1.append(temp_str)
    # for key,val in str1:
    #     res_matrix.setdefault(key,[]).append(val)


    # calculate polymorphics sites format
    nuclOrder = ['A', 'C', 'G', 'T']
    polymorphicSites = []
    for loc in range(len(res)):
        pos_nucl = res[loc]
        formatOutput = [loc]

        flag_nonZero = 0
        for nucl in nuclOrder:
            nucl_count = pos_nucl[nucl]
            formatOutput.append(nucl_count)
            if nucl_count != 0:
                flag_nonZero += 1

        if flag_nonZero >= 1:
            polymorphicSites.append(formatOutput)



    # res_file_name = res_dir + '/merged_polymorphic_sites_without_filter'
    # with open('%s' % res_file_name, 'w') as f:
    #     for item in polymorphicSites:
    #         one = '\t'.join([str(one) for one in item])
    #         f.write('%s\n' % one)

    # a nucleotide at each position is from fewer than 5% of samples or no larger than 5 or fewer than 5% of the total coverage of that position
    polymorphicSites_filter1 = []
    for item in polymorphicSites:
        temp_item = [item[0]]
        for one in item[1:]:
            index1 = item[1:].index(one)
            sample_num = len(res_count[item[0]][nuclOrder[index1]])
            if one <= 5 or sample_num <= len(reads_file_list) * 0.05 or one <= sum(item[1:])*0.05:
                temp_item.append(0)
            else:
                temp_item.append(one)
        polymorphicSites_filter1.append(temp_item)

    # res_file_name1 = res_dir + '/merged_polymorphic_sites_filter_2B'
    # with open('%s' % res_file_name1, 'w') as f:
    #     for item in polymorphicSites_filter1:
    #         one = '\t'.join([str(one) for one in item])
    #         f.write('%s\n' % one)

    ########## updated reference genome ###############
    updated_genome(polymorphicSites_filter1, genomeFileLoc)

    ######### filter four zero ########################
    polymorphicSites_filter2 = filter_four_zero(polymorphicSites_filter1)

    ###### filter polymorphic sites less than 10% of the average coverage ####
    polymorphicSites_filter3 = filter(polymorphicSites_filter2)

    res_file_name = res_dir + '/merged_filter_polymorphic_sites'
    with open('%s' % res_file_name, 'w') as f:
        for item in polymorphicSites_filter3:
            one = '\t'.join([str(one) for one in item])
            f.write('%s\n' % one)

    # res_file_name1 = res_dir + '/input_each_sample_polymorphic_sites'
    # with open('%s' % res_file_name1, 'w') as f:
    #     for key,val in res_matrix.items():
    #         f.write(str(key) + '\t')
    #         one = '\t'.join([str(one) for one in val])
    #         f.write('%s\n' % one)
    # return polymorphicSites_filter3,res_matrix
    return polymorphicSites_filter3

# genome_len = 2821361
# res_tmp_dir = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/Min_data/chang_paper/Aureus/NC_007795.1'
# reads_file_list = reads_frequency_file_list = [res_tmp_dir + '/sample_reads_frequency/' + i for i in os.listdir(res_tmp_dir + '/sample_reads_frequency')[:3]]
# res_dir ='/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/Min_data/chang_paper/Aureus'
# genomeFileLoc = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/Min_data/chang_paper/original_genome/AureusIndex/GCF_000013425.1_ASM1342v1_genomic_RefAureus.fna'
# getMergedPolymorphicSites(genome_len,reads_file_list, res_dir,genomeFileLoc)