import os
from copy import deepcopy

REF_COUNT_NORMALIZED = 5
SECOND_COUNT_NORMALIZED = 5


def read_polymorphicsSite(inName):
    data = []
    with open('%s' % inName) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            line = [int(one) for one in line]
            data.append(line)
    return data

def check_refernce_genome(polyCounts, genomeFileLoc):
    genome = ''
    with open('%s' % genomeFileLoc) as f:
        for line in f:
            if line.startswith('>'):
                pass
            else:
                genome += line.strip()

    polyCounts_rem = {}
    order = ['A', 'C', 'G', 'T']
    for item in polyCounts:
        pos = item[0]
        ref = genome[pos].upper()
        polyCounts_rem.update({pos: ref})

    return polyCounts_rem

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

# def get_max_sum_frequency(data_raw):
#
#     list_all = []
#     for i in range(len(data_raw)):
#         tmp = sum(data_raw[i][1:5])
#         list_all.append(tmp)
#
#     max_value = max(list_all)
#
#     if max_value > 100:
#         return max_value
#     else:
#         return 100
import numpy as np
# def get_median_sum_frequency(data_raw):
#
#     list_all = []
#     for i in range(len(data_raw)):
#         tmp = sum(data_raw[i][1:5])
#         list_all.append(tmp)
#
#     max_value = np.median(list_all)
#
#     # if max_value > 100:
#     #     return max_value
#     # else:
#     #     return 100
#
#     return max_value


def normalized_to_max(count,max_value):
    freq = []
    coe = max_value / sum(count)
    for i in range(len(count)):
        freq.append(int(round(count[i] * coe, 0)))

    # make it is A+C+G+T=100, may be not necessay
    if sum(freq) != max_value:
        difference = max_value - sum(freq)
        maxIndex = [index for index, value in enumerate(freq) if
                    value == max(freq)]
        freq[maxIndex[0]] += difference

    return freq

def normalized_filter(data, data_ref):
    data_raw = deepcopy(data)
    data_raw = filter_four_zero(data_raw)
    # normalized polymorphic sites to sum of 100. A+C+G+T=100. and only pick
    # the data in the format = [ref, ref_index]
    polyCounts = []
    orderMap = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # max_value = get_max_sum_frequency(data_raw)
    max_value = 100
    for i in range(len(data_raw)):
        tmp = normalized_to_max(data_raw[i][1:5],max_value)

        # get ref
        loc = data_raw[i][0]
        ref_index = orderMap[data_ref[loc]]

        ref_count = tmp[ref_index]
        # filter ref count
        # if ref_count < REF_COUNT_NORMALIZED:
        #     continue
        # the largest number
        # first_count = sorted(tmp)[-1]
        # if first_count <= SECOND_COUNT_NORMALIZED:
        #     continue
        # the second largest number.
        second_count = sorted(tmp)[-2]

        if second_count <= SECOND_COUNT_NORMALIZED:
            continue
        tmp = [data_raw[i][0], ref_count,data_ref[loc]] + tmp

        polyCounts.append(tmp[:])
    polyCounts.sort(key=lambda item: (item[1], item[2]), reverse=True)
    return polyCounts

def read_polymorphicsSite_sub(inName):
    data = []
    with open('%s' % inName) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            data.append(line)
    return data

def normalized_matrix(inFileName, genomeFileLoc):
    # read polymorphic sites and filter
    data = read_polymorphicsSite(inFileName)

    # label polymorphic sites what is reference in that pos.
    polyCounts_refLabel = check_refernce_genome(data, genomeFileLoc)

    # normalized and filter individual letter
    print('normalize and filter data')
    data_normalized = normalized_filter(data, polyCounts_refLabel)

    loc_list = []
    with open('%s' % inFileName + '_normalized_100', 'w') as f:
        for line in data_normalized:
            loc_list.append(line[0])
            tmp = '\t'.join([str(one) for one in line])

            # tmp = '\t'.join([str(line[0]),str(line[3]),str(line[4]),str(line[5]),str(line[6])])

            f.write('%s\n' % tmp)
    # resName = os.path.split(inFileName)[0] + '/input_each_sample_polymorphic_sites'
    resName1 = inFileName + '_100'
    # data1_matrix = read_polymorphicsSite_sub(resName)
    # print(loc_list)
    # print(data1_matrix)
    with open('%s' % resName1, 'w') as f:
        for line in data:
            if line[0] in loc_list:

                fre = sorted(line[1:])

                line1 = [line[0]]

                for i in line[1:]:
                    if i < fre[-2]:
                        line1.append(0)
                    else:
                        line1.append(i)

                tmp = '\t'.join([str(one) for one in line1])
                f.write('%s\n' % tmp)