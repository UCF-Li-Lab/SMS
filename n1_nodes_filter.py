
import os


def get_nodes(merged_file):
    f1 = open(merged_file, 'r')
    lines1 = f1.readlines()
    nodes_list = []
    ## nucleotid_position_frequency
    # order = ['A', 'C', 'G', 'T']
    for line1 in lines1:
        loc = line1.split('\t')[0]
        ref = line1.split('\t')[2]
        A_fre = line1.split('\t')[3]
        C_fre = line1.split('\t')[4]
        G_fre = line1.split('\t')[5]
        T_fre = line1.split('\t')[6].strip()

        if A_fre != '0' and ref != 'A':
            node = 'A' + '_' + loc + '_' + A_fre
            nodes_list.append(node)

        if C_fre != '0' and ref != 'C':
            node = 'C' + '_' + loc + '_' + C_fre
            nodes_list.append(node)

        if G_fre != '0' and ref != 'G':
            node = 'G' + '_' + loc + '_' + G_fre
            nodes_list.append(node)

        if T_fre != '0' and ref != 'T':
            node = 'T' + '_' + loc + '_' + T_fre
            nodes_list.append(node)

    return nodes_list

def read_polymorphicsSite(inName):
    data = {}
    with open('%s' % inName) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            data[line[0]] = [line[1], line[2], line[3], line[4]]
    return data

def normalized_to_max(count, max_value):
    freq = []
    coe = max_value / sum(count)
    for i in range(len(count)):
        freq.append(round(count[i] * coe))
    return freq

def get_nodes_frequency_vectors(n1_nodes, frequency_file_loc):
    reads_file_list = [frequency_file_loc + '/' + i for i in os.listdir(frequency_file_loc)]

    print('generaring node vector, this will take some time')

    dict1 = {}
    str1 = []
    nodes_list = n1_nodes
    node_loc_list = list(set([i.split('_')[1] for i in nodes_list]))
    for reads_file in reads_file_list:
        # generate polymorphic sites
        print('Reading file---------' + str(reads_file_list.index(reads_file)))
        # print(reads_file)
        sub_polymorphicSites = read_polymorphicsSite(reads_file)

        str1 += [[node, sub_polymorphicSites[node]] for node in node_loc_list]

    for key, val in str1:
        dict1.setdefault(key, []).append(val)

    dict2 = {}  #### normalized frequency vector dict
    dict2_1 = {}  #### original frequency vector dict
    order_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    normalized_dict = {}
    print(n1_nodes)
    for node in n1_nodes:
        L = node.split('_')[0]
        loc = node.split('_')[1]

        normalized_fre = int(node.split('_')[2])

        normalized_dict[L + '_' + loc] = normalized_fre

        vector = dict1[loc]
        vector1 = []
        for i in vector:
            vector1.append(int(i[order_dict[L]]))

        normalized_vector1 = normalized_to_max(vector1, normalized_fre)

        if sum(normalized_vector1) >=6:
            dict2[node] = normalized_vector1
    return dict2


def get_file(loc_list,file1,file2):
    f1 = open(file1,'r')
    f2 = open(file2,'w')

    lines1 = f1.readlines()

    for line1 in lines1:

        loc = line1.split('\t')[0]

        if loc in loc_list:
            f2.writelines(line1)


    f1.close()
    f2.close()



def get_n1_nodes_filter(merge_file,frequency_file_loc,result_dir):
    print('getting vectors of the nodes, this will cost some time')

    nodes_list = get_nodes(merge_file)
    node_loc_list = list(set([str(int(i.split('_')[1])+1) for i in nodes_list]))
    #
    with open(result_dir + '/all_nodes_loc','w') as f:
        f.write('\n'.join(node_loc_list))
    # # # # # # #
    reads_file_list = os.listdir(frequency_file_loc)
    # # # # # # # #
    filter_reads_file = frequency_file_loc + '_filter'
    # # # # # # # #
    if os.path.exists(filter_reads_file) == False:
        os.makedirs(filter_reads_file)
    #
    for file in reads_file_list:
    #
        print('Filter file -----------',reads_file_list.index(file))
    #
        input_file = frequency_file_loc + '/'+ file
        output_file = filter_reads_file + '/' + file
        node_file = result_dir + '/all_nodes_loc'
    #
        command1 = "awk 'FNR==NR { seen[$0]++ }; FNR!=NR && FNR in seen' %s %s > %s" % (node_file,input_file,output_file)
    #
        os.system(command1)

    dict_all = get_nodes_frequency_vectors(nodes_list, filter_reads_file)

    input_file1 = merge_file
    output_file1 = merge_file + '_filter'

    input_file2 = result_dir + '/merged_filter_polymorphic_sites_100'
    output_file2 = result_dir + '/merged_filter_polymorphic_sites_100_filter'

    loc_list = [i.split('_')[1] for i in list(dict_all.keys())]

    get_file(loc_list,input_file1,output_file1)
    get_file(loc_list,input_file2,output_file2)


# merge_file = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/xin_195_data/final_result_on_min_website/NC_000962.3/merged_filter_polymorphic_sites_normalized_100'
# frequency_file_loc = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/xin_195_data/final_result_on_min_website/NC_000962.3/sample_reads_frequency'
# result_dir = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/xin_195_data/final_result_on_min_website/NC_000962.3'
#
# get_n1_nodes_filter(merge_file,frequency_file_loc,result_dir)