import os
import warnings
from random import random

from scipy.stats import poisson
import math
from copy import deepcopy
import numpy as np
import Bio.Cluster
from collections import Counter


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


def normalized_to_max(count, max_value):
    freq = []
    coe = max_value / sum(count)
    for i in range(len(count)):
        freq.append(round(count[i] * coe))
    return freq

def normalized_to_max1(count, max_value):
    freq = []
    coe = max_value / sum(count)
    for i in range(len(count)):
        freq.append(round(count[i] * coe,2))
    return freq


def read_polymorphicsSite(inName):
    data = {}
    with open('%s' % inName) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            data[line[0]] = [line[1], line[2], line[3], line[4]]
    return data


def get_nodes_frequency_vectors(n1_nodes, R_nodes, frequency_file_loc):
    reads_file_list = [frequency_file_loc + '/' + i for i in os.listdir(frequency_file_loc)]

    print('generaring node vector, this will take some time')

    dict1 = {}
    str1 = []
    nodes_list = n1_nodes + R_nodes
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

        dict2[node] = normalized_vector1
        dict2_1[node] = vector1

    dict3 = {}
    order_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for node in R_nodes:
        L = node.split('_')[0]
        loc = node.split('_')[1]

        vector = dict1[loc]
        vector1 = []
        for i in vector:
            vector1.append(int(i[order_dict[L]]))

        if L + '_' + loc in normalized_dict.keys():
            normalized_fre = normalized_dict[L + '_' + loc]

            normalized_vector1 = normalized_to_max(vector1, normalized_fre)
            dict3[node] = normalized_vector1

    return [dict2, dict3, dict2_1]


def get_kendalltau_correlaton(vector1, vector2):
    kendalltau_score = 1 - Bio.Cluster.distancematrix((vector1, vector2), dist="k")[1][0]
    return kendalltau_score


def poisson_pmf(x, mu):
    try:
        return mu ** x / math.factorial(x) * math.exp(-mu)
    except:
        return stats.poisson.pmf(x,mu)


#### difference of the sum of the two vectores ######

def vector_sum_diff(vector1, vector2):
    a = round(sum([i for i in vector1]))
    b = round(sum([i for i in vector2]))
    value_abs = poisson_pmf(a, b) + poisson_pmf(b, a)
    # value_abs = stats.poisson.pmf(a,b) + stats.poisson.pmf(b,a)
    return value_abs


def prediction(sum_diff, kendalltau):
    score = 79.25 * sum_diff + 43.06 * kendalltau + 43.06 * math.pow(kendalltau, 3) + 0.04 / (
                -0.0025 - sum_diff) - 11.45
    return score


def similar_formula(vector1, vector2):
    sum_diff = vector_sum_diff(vector1, vector2)
    kendalltau_score = get_kendalltau_correlaton(vector1, vector2)
    score = prediction(sum_diff, kendalltau_score)
    return score


def get_R_strain(file1):
    f1 = open(file1, 'r')
    lines1 = f1.readlines()
    lines1 = ''.join(lines1).split('>')
    dict1 = {}

    for line1 in lines1[1:]:
        index1 = lines1.index(line1)
        cov = line1.split('\n')[0]
        node_list = [i.split(',')[0] + '_' + i.split(',')[1] for i in line1.split('\n')[1:-1]]
        dict1[str(index1) + '_' + cov] = node_list

    return dict1


def get_strain_average_vector(nodes_list, nodes_vector_dict, vector_size):
    vector_average = []
    for i in range(vector_size):
        temp_list = []
        for node in nodes_list:
            if node in nodes_vector_dict.keys():
                temp_list.append(int(nodes_vector_dict[node][i]))

        ave = float(sum(temp_list)) / float(len(temp_list))
        ave1 = round(float(sum(temp_list)) / float(len(temp_list)))

        ######### in case the averge is too small #################
        # if vector_size > 50:
        #     vector_average.append(ave)
        # else:
        vector_average.append(ave1)
        # vector_average.append(ave)
    # print('vector average')
    # print(vector_average)

    return vector_average


def get_R_strain_vector(R_strain_dict, nodes_vector_dict):
    R_average_vector = {}
    for key, val in R_strain_dict.items():
        nodes_list = val
        vector_size = len(list(nodes_vector_dict.values())[0])
        strain_average_vector = get_strain_average_vector(nodes_list, nodes_vector_dict, vector_size)
        R_average_vector[key] = strain_average_vector
    return R_average_vector


def get_new_strain_index(n1_nodes_vector_dict, strain_vector):
    score_list = []
    for key, val in n1_nodes_vector_dict.items():
        score_temp = []
        y = val
        for key1, val1 in strain_vector.items():
            x = val1
            score = similar_formula(x, y)
            score_temp.append(score)
        score_list.append(score_temp)

    n1_new_strain_index = []
    for score in score_list:
        score_index = sorted(range(len(score)), key=lambda k: score[k])
        list1 = [score_index[-1]]  ######## max score
        n1_new_strain_index.append(list(set(list1)))
        # print(score,list(set(list1)))
    return n1_new_strain_index


def get_n1_nodes_strain(n1_nodes_vector_dict, R_strain_vector_dict, file1):
    n1_strain_vector_dict = {}
    str1 = []
    n1_new_strain_index = get_new_strain_index(n1_nodes_vector_dict, R_strain_vector_dict)
    for i in range(len(n1_new_strain_index)):
        index1_list = n1_new_strain_index[i]
        for index1 in index1_list:
            new_key = list(R_strain_vector_dict.keys())[index1]
            str1.append([new_key, list(n1_nodes_vector_dict.keys())[i]])

    for key, val in str1:
        n1_strain_vector_dict.setdefault(key, []).append(val)

    previous_n1_nodes_strain = ''
    ROUND = 1

    f1 = open(file1, 'w')
    while True:
        strain_vector = get_R_strain_vector(n1_strain_vector_dict, n1_nodes_vector_dict)
        str1 = []
        n1_new_strain_index = get_new_strain_index(n1_nodes_vector_dict, strain_vector)
        for i in range(len(n1_new_strain_index)):
            index1_list = n1_new_strain_index[i]
            for index1 in index1_list:
                new_key = list(R_strain_vector_dict.keys())[index1]
                str1.append([new_key, list(n1_nodes_vector_dict.keys())[i]])

        n1_strain_vector_dict = {}
        for key, val in str1:
            n1_strain_vector_dict.setdefault(key, []).append(val)

        # check when the ROUND is stopped
        flag_converge = True
        if ROUND == 1:
            flag_converge = False
            previous_n1_nodes_strain = deepcopy(n1_strain_vector_dict)
        else:

            change_num = 0
            for key, val in n1_strain_vector_dict.items():
                # print(key,previous_n1_nodes_strain.keys())
                if key in previous_n1_nodes_strain.keys():
                    num1 = len(list(set(val).symmetric_difference(set(previous_n1_nodes_strain[key]))))
                    change_num += num1
            print(change_num)
            if change_num > 0:
                flag_converge = False

            previous_n1_nodes_strain = deepcopy(n1_strain_vector_dict)

            f1.write(str(ROUND) + '\t' + str(change_num) + '\n')

        if ROUND > 200:
            flag_converge = True
            f1.close()

        if flag_converge:

            n1_strain_vector_dict1 = {}

            for key, val in n1_strain_vector_dict.items():
                f_val = [int(i.split('_')[-1]) for i in list(set(val))]
                ave_val = round(float(sum(f_val)) / float(len(f_val)), 2)
                new_key = key.split('_')[0] + '_' + str(ave_val)
                n1_strain_vector_dict1[new_key] = list(set(val))
            f1.close()
            return n1_strain_vector_dict1

        ROUND += 1

        print('ROUND----------------------------' + str(ROUND))


def get_nodes_assigned(val_dict, R_strain_vector_dict, current_strain_dict):
    for key1, val1 in val_dict.items():
        dis_list = []
        for key2, val2 in R_strain_vector_dict.items():
            dis = similar_formula(val1, val2)
            dis_list.append(dis)

        index_max = dis_list.index(max(dis_list))
        new_key = list(R_strain_vector_dict.keys())[index_max]
        current_strain_dict[new_key].append(key1)

    return current_strain_dict


from scipy import stats


def ZIP(x, pi, mu):
    return (x == 0) * pi + (1 - pi) * stats.poisson.pmf(x, mu)


from sklearn.mixture import GaussianMixture
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning,
                        module="sklearn")


# def divide_into_two_group(node_dict):
#     X = np.array(list(node_dict.values()))
#
#     # gm = GaussianMixture(n_components=2,random_state=10).fit(X)
#     # print('X')
#     # print(node_dict.keys())
#     gm = GaussianMixture(n_components=2,random_state=10,n_init=20).fit(X)
#     labels = list(gm.predict(X))
#     dict1 = {}
#     str1 = []
#     for i in range(len(labels)):
#         label = str(labels[i])
#         node = list(node_dict.keys())[i]
#         str1.append([label, node])
#
#     for key, val in str1:
#         dict1.setdefault(key, []).append(val)
#
#     group1 = list(dict1.values())[0]
#     group2 = list(dict1.values())[1]
#     print('group1')
#     # print(group1)
#     print('group2')
#     # print(group2)
#     return [group1, group2]

def divide_into_two_group(node_dict):
    X = np.array(list(node_dict.values()))

    # gm = GaussianMixture(n_components=2,random_state=10).fit(X)
    # print('X')
    # print(node_dict.keys())
    gm = GaussianMixture(n_components=2, random_state=10, n_init=20).fit(X)
    labels = list(gm.predict(X))
    dict1 = {}
    str1 = []
    for i in range(len(labels)):
        label = str(labels[i])
        node = list(node_dict.keys())[i]
        str1.append([label, node])

    # print('labels')
    # print(labels)

    for key, val in str1:
        dict1.setdefault(key, []).append(val)

    if len(list(dict1.keys())) == 1:
        group1 = list(dict1.values())[0]
        group2 = [group1[0]]
        return [group1, group2]

    else:

        group1 = list(dict1.values())[0]
        group2 = list(dict1.values())[1]
        print('group1')
        # print(group1)
        print('group2')
        # print(group2)
        return [group1, group2]


# def divide_into_two_group(node_dict):
#     X = np.array(list(node_dict.values()))
#
#     # aic1 = GaussianMixture(n_components=1, random_state=0).fit(X).aic(X)
#     # bic1 = GaussianMixture(n_components=1, random_state=0).fit(X).bic(X)
#
#     # gm = GaussianMixture(n_components=2,n_init=2 ,random_state=0,covariance_type='tied').fit(X)
#     gm = GaussianMixture(n_components=2, random_state=0).fit(X)
#
#     # aic = GaussianMixture(n_components=2, random_state=0).fit(X).aic(X)
#     # bic = GaussianMixture(n_components=2, random_state=0).fit(X).bic(X)
#
#     gm = GaussianMixture(n_components=2, random_state=10).fit(X)
#     labels = list(gm.predict(X))
#     dict1 = {}
#     str1 = []
#     for i in range(len(labels)):
#         label = str(labels[i])
#         node = list(node_dict.keys())[i]
#         str1.append([label, node])
#
#     for key, val in str1:
#         dict1.setdefault(key, []).append(val)
#
#     group1 = list(dict1.values())[0]
#     group2 = list(dict1.values())[1]
#     print('group1')
#     print(group1)
#     print('group2')
#     print(group2)
#     return [group1, group2]


import scipy.stats


def BIC(likelihood_log, k, data_point_num):
    bic = k * math.log(data_point_num) - 2 * likelihood_log
    return bic


def get_score1(likehood_diff,deg_diff):
    X = 2 * likehood_diff
    k = deg_diff

    if k != 0:
        x = (X-k)/math.pow(2*k,0.5)
        return x
    else:
        return 1000 ########### we set cutoff 300 to divide, if k ==0, we will divide it ######

def get_score2(likehood_diff,deg_diff):
    X = 2 * likehood_diff
    k = deg_diff
    x = math.pow(2*X,0.5)
    # mean = math.pow(2*k-1,0.5)

    # pdf = norm.pdf(x,loc=mean,scale=1)

    return x


def get_score3(likehood_diff, deg_diff):
    X = 2 * likehood_diff
    k = deg_diff

    if k != 0:
        x = math.pow(X/k, 1/3)
    # mean = 1- 2/(9*k)
    # scale1 = math.pow(2/(9*k),0.5)
    # pdf = norm.pdf(x, loc=mean, scale=scale1)

        return x
    else:
        return 1000


def get_likelihood_ratio(diff, d):
    # print(diff,d)
    p_value = scipy.stats.chi2.sf(2 * diff, d)
    return p_value


def covert_na(a1):
    # np.place(a1, a1 == 'na', np.nan)
    # print('using np.where')
    # print(a1)
    # a1[a1 =='na'] = np.nan
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=FutureWarning)
        a1 = np.where(a1 == 'na', np.nan, a1)
        a1 = a1.astype(np.float)

        return a1


def Newton(b_list, n_list, a_list, x0_list, strain_num):
    from scipy import optimize

    # print(b_list)
    # print(n_list)
    # print(a_list)
    # print(x0_list)
    # print(strain_num)

    root_list = []
    for r in range(strain_num):
        def f(x):
            with np.errstate(all='ignore'):
                pi_r_val1 = n_list[:, r] * (x - a_list[:, r]) * np.exp(-x)
                pr_r_val2 = x * b_list[:, r] - n_list[:, r] * (x - a_list[:, r]) * (1 - np.exp(-x))
                pi_r = pi_r_val1 / pr_r_val2

                f_val1 = n_list[:, r] * a_list[:, r] / x
                f_val2 = (1 - pi_r) * b_list[:, r] / (pi_r + (1 - pi_r) * np.exp(-x))

                f = f_val1 - f_val2

            return f

        with np.errstate(all='ignore'):
            root = optimize.newton(f, x0_list[:, r], maxiter=1000)
            root_list.append(root)
    return np.array(root_list).T


################### step8 ######################
def get_likehood(strain_dict, nodes_vector_dict):
    strain_num = len(strain_dict.keys())
    sample_num = len(list(nodes_vector_dict.values())[0])

    def get_X_list(nodes_vecotr_dict):
        X_list = []

        for key, val in nodes_vecotr_dict.items():
            X_list.append(val)

        return X_list

    def convert_T(a, b):
        if a == b:
            return 1
        else:
            return 0

    X_list = np.array(get_X_list(nodes_vector_dict))
    covert = np.vectorize(convert_T)
    X_list_T = covert(X_list, 0)

    def get_z_list(strain_dict, nodes_vector_dict):
        z_list = []
        for i in range(len(nodes_vector_dict.keys())):
            z_temp = []
            node1 = list(nodes_vector_dict.keys())[i]
            for r in range(strain_num):
                r_nodes_list = set(list(strain_dict.values())[r])
                if node1 in r_nodes_list:
                    z_temp.append(1)
                else:
                    z_temp.append(0)

            z_list.append(z_temp)

        return z_list

    z_list = np.array(get_z_list(strain_dict, nodes_vector_dict))

    def get_Y_list(strain_dict, nodes_vector_dict):
        Y_list = []
        Y_list_non_zero = []
        for j in range(sample_num):
            y_temp = []
            y_temp1 = []
            for r in range(strain_num):
                r_nodes_list = set(list(strain_dict.values())[r])
                non_zero_list = []
                y_temp1_temp = 0
                for r_nodes_list_temp in r_nodes_list:
                    if nodes_vector_dict[r_nodes_list_temp][j] > 0:
                        non_zero_list.append(1)
                        y_temp1_temp += 1

                if len(non_zero_list) >= 5:
                    y_temp.append(1)

                else:
                    y_temp.append(0)

                y_temp1.append(y_temp1_temp)
            Y_list.append(y_temp)
            Y_list_non_zero.append(y_temp1)

        return [Y_list, Y_list_non_zero]

    Y_list_all = get_Y_list(strain_dict, nodes_vector_dict)

    Y_list = np.array(Y_list_all[0])
    Y_list_non_zero = Y_list_all[1]

    def get_a_b_n(z_list, X_list, X_list_T, Y_list, strain_num):

        a_list = []
        b_list = []
        n_list = []

        for j in range(sample_num):
            b_temp = []
            n_temp = []
            a_temp = []
            for r in range(strain_num):
                b_num = sum(X_list_T[:, j] * z_list[:, r]) * Y_list[j][r]
                n_num = sum(z_list[:, r] * Y_list[j][r])

                if n_num != 0:
                    n_temp.append(n_num)
                    a_num = sum(z_list[:, r] * X_list[:, j]) * Y_list[j][r] / n_num
                    if a_num == 0:
                        a_temp.append('na')
                    else:
                        a_temp.append(round(a_num, 3))
                else:
                    n_temp.append('na')
                    a_num = 'na'
                    a_temp.append(a_num)

                if b_num == 0:
                    b_temp.append('na')
                    # b_temp.append(1e-05)
                    # b_temp.append(1)
                else:
                    b_temp.append(b_num)

            a_list.append(a_temp)
            b_list.append(b_temp)
            n_list.append(n_temp)
        return [a_list, b_list, n_list]

    a_b_n_list = get_a_b_n(z_list, X_list, X_list_T, Y_list, strain_num)
    a_list = a_b_n_list[0]
    b_list = a_b_n_list[1]
    n_list = a_b_n_list[2]

    def get_s_pow_initional(X_list, a_list, z_list):
        s_pow_initional_list = []
        for j in range(sample_num):
            s_temp = []
            for r in range(strain_num):
                if a_list[j][r] != 'na':
                    val1 = sum(z_list[:, r] * Y_list[j][r] * np.square(X_list[:, j] - a_list[j][r]))
                    val2 = sum(z_list[:, r])

                    if val2 == 1:
                        s_temp.append('na')
                    else:
                        s_r = val1 / (val2 - 1)
                        s_temp.append(round(s_r, 3))
                else:
                    s_temp.append('na')
            s_pow_initional_list.append(s_temp)

        return s_pow_initional_list

    s_pow_initional_list = get_s_pow_initional(X_list, a_list, z_list)

    def get_lamda_initional(s_pow_initional_list):

        lamda_initional_list = []
        for j in range(sample_num):
            lamda_temp = []
            for r in range(strain_num):
                if s_pow_initional_list[j][r] != 'na' and a_list[j][r] != 'na':
                    val1 = s_pow_initional_list[j][r] + math.pow(a_list[j][r], 2)
                    val2 = a_list[j][r]
                    with np.errstate(all='ignore'):
                        lamda_val = val1 / val2 - 1
                        lamda_temp.append(round(lamda_val, 3))
                else:
                    lamda_temp.append('na')

            lamda_initional_list.append(lamda_temp)

        return lamda_initional_list

    lamda_initional_list = get_lamda_initional(s_pow_initional_list)

    def get_lamda_pi(lamda_initional_list, a_list, n_list, b_list, strain_num):

        b_list = covert_na(b_list)
        n_list = covert_na(n_list)
        a_list = covert_na(a_list)
        lamda_initional_list = covert_na(lamda_initional_list)

        lamda_list = Newton(b_list, n_list, a_list, lamda_initional_list, strain_num)

        pi_list = []
        for r in range(strain_num):
            pi_r_val1 = n_list[:, r] * (lamda_list[:, r] - a_list[:, r]) * np.exp(-lamda_list[:, r])
            pr_r_val2 = lamda_list[:, r] * b_list[:, r] - n_list[:, r] * (lamda_list[:, r] - a_list[:, r]) * (
                        1 - np.exp(-lamda_list[:, r]))
            pi_r = pi_r_val1 / pr_r_val2

            pi_list.append(pi_r)

        pi_list = np.array(pi_list).T

        return [lamda_list, pi_list]

    lamda_pi_list = get_lamda_pi(np.array(lamda_initional_list), np.array(a_list), np.array(n_list), np.array(b_list),
                                 strain_num)
    lamda_list = lamda_pi_list[0]
    pi_list = lamda_pi_list[1]

    ####### convert to zero if pi<1 or pi >1
    pi_list = np.where(pi_list < 0, 0, pi_list)
    pi_list = np.where(pi_list > 1, 0, pi_list)
    pi_list = pi_list.astype(np.float)

    mask = (pi_list == 0)
    new_array = np.copy(lamda_list)

    # print(a_list)

    new_array[mask] = np.array(a_list)[mask]
    lamda_list = new_array

    # print('lamda_list')
    # print(lamda_list)
    # print('pi_list')
    # print(pi_list)

    # for i in range(sample_num):
    #     for r in range(strain_num):
    #         if Y_list[i][r] == 1 and np.isnan(pi_list[i][r])== True:
    #             pi_list[i][r] = 0
    #             lamda_list[i][r] = a_list[i][r]

    mask1 = (np.isnan(pi_list))
    mask2 = (Y_list == 1)
    mask3 = mask1 & mask2
    # new_array1 = np.copy(lamda_list)
    lamda_list[mask3] = np.array(a_list)[mask3]
    # lamda_list = new_array1
    pi_list[mask3] = 0

    def log_likelihood(X_list, pi_list, lamda_list, z_list, Y_list):
        log_likehood = 0
        current_log_list = []
        node_num = np.shape(X_list)[0]

        for i in range(node_num):
            temp1 = []
            for r in range(strain_num):
                x_i = X_list[i, :]
                pi_r = pi_list[:, r]
                lamda_r = lamda_list[:, r]

                if z_list[i][r] == 0:
                    val1 = np.zeros(sample_num)
                    temp1.append(val1)
                else:
                    val = ZIP(x_i, pi_r, lamda_r)
                    val1 = z_list[i][r] * Y_list[:, r] * val
                    val1 = np.nan_to_num(val1, nan=0)
                    temp1.append(val1)

            temp2 = sum(temp1)
            log_sum = sum(np.log([i for i in temp2 if i > 0]))
            current_log_list.append(log_sum)
            log_likehood += log_sum

        return log_likehood

    log_likehood = log_likelihood(X_list, pi_list, lamda_list, z_list, Y_list)

    return [log_likehood, pi_list, Y_list_non_zero]


def save_data(data, file1):
    import pickle
    file = open(file1, 'wb')

    # dump information to that file
    pickle.dump(data, file)

    # close the file
    file.close()


def load_data(file1):
    import pickle

    # open a file, where you stored the pickled data
    file = open(file1, 'rb')

    # dump information to that file
    data = pickle.load(file)

    # close the file
    file.close()

    return data


def get_parameter_number(pi_list):
    a = np.count_nonzero(pi_list > 0)
    b = np.count_nonzero(pi_list == 0)
    num1 = 2 * a + b

    return num1


########### R and R-1 step data point ###########
def get_data_point_number(Y_non_zero_list):
    num1 = 0
    for i in Y_non_zero_list:
        num1 += sum(i)
    return num1


########## R and R+1 step data point #############
def get_data_point_number1(Y_non_zero_list, nodes_number_list):
    num1_list = []

    for j in range(len(Y_non_zero_list[0])):
        temp_num = []
        for i in Y_non_zero_list:
            if i[j] != 0:
                temp_num.append(i[j])

        num1_list.append(temp_num)

    num_all = 0
    # print(num1_list)

    for i in range(len(num1_list)):
        num_all += nodes_number_list[i] * len(num1_list[i])
    return num_all


########## check R and R-1 strain ##################################
def check_R_strain(n1_nodes_strain_dict, n1_nodes_vector_dict, BIC_file):
    previous_strain_dict = deepcopy(n1_nodes_strain_dict)
    # for key in sorted(n1_nodes_strain_dict.keys(), key=lambda x: len(n1_nodes_strain_dict[x])):
    f1 = open(BIC_file, 'w')
    key_list = sorted(previous_strain_dict.keys(), key=lambda x: len(previous_strain_dict[x]))

    key = key_list[0]
    val = previous_strain_dict[key]
    ROUND = 0
    while True:
        previous_result = get_likehood(previous_strain_dict, n1_nodes_vector_dict)
        previous_likehood = previous_result[0]
        previous_pi_list = previous_result[1]
        previous_non_zero_list = previous_result[2]

        current_strain_dict = deepcopy(previous_strain_dict)
        del current_strain_dict[key]

        val_dict = {}

        for i in val:
            val_dict[i] = n1_nodes_vector_dict[i]

        strain_vector_dict = get_R_strain_vector(current_strain_dict, n1_nodes_vector_dict)
        current_strain_dict = get_nodes_assigned(val_dict, strain_vector_dict, current_strain_dict)

        current_result = get_likehood(current_strain_dict, n1_nodes_vector_dict)
        current_likehood = current_result[0]
        current_pi_list = current_result[1]
        current_non_zero_list = current_result[2]

        previous_parameter_num = get_parameter_number(previous_pi_list)

        previous_point_num = get_data_point_number(previous_non_zero_list)

        previous_bic = BIC(previous_likehood, previous_parameter_num, previous_point_num)

        current_parameter_num = get_parameter_number(current_pi_list)

        current_point_num = get_data_point_number(current_non_zero_list)

        current_bic = BIC(current_likehood, current_parameter_num, current_point_num)

        diff = abs(previous_likehood - current_likehood)
        deg_diff = abs(previous_parameter_num - current_parameter_num)

        likelihood_p = get_likelihood_ratio(diff, deg_diff)

        f1.write(str(previous_strain_dict.keys()) + '\t' + str(current_strain_dict.keys()) + '\t' + str(
            previous_bic) + '\t' + str(current_bic) + '\t' + str(
            [previous_likehood, previous_parameter_num, previous_point_num]) + '\t'
                 + str([current_likehood, current_parameter_num, current_point_num]) + '\n')

        score1 = get_score1(diff,deg_diff)
        score2 = get_score2(diff, deg_diff)
        score3 = get_score3(diff, deg_diff)

        print('likelihood_p')
        print(likelihood_p, deg_diff, diff, score1,score2,score3)

        # bic_diff = previous_bic - current_bic
        # label = False
        # if deg_diff == 0 and bic_diff > 0:
        #     label = True

        if likelihood_p > 1e-110:
            previous_strain_dict = current_strain_dict
            key_list = sorted(previous_strain_dict.keys(), key=lambda x: len(previous_strain_dict[x]))
            key = key_list[0]
            val = previous_strain_dict[key]

            if len(previous_strain_dict) == 1:
                f1.close()
                previous_strain_dict1 = {}

                for key, val in previous_strain_dict.items():
                    f_val = [int(i.split('_')[-1]) for i in val]
                    ave_val = round(float(sum(f_val)) / float(len(f_val)), 2)
                    index1 = list(previous_strain_dict.keys()).index(key)
                    new_key = str(index1) + '_' + str(ave_val)
                    previous_strain_dict1[new_key] = val

                return previous_strain_dict1
        else:

            if ROUND >= len(key_list) - 1:
                f1.close()

                previous_strain_dict1 = {}

                for key, val in previous_strain_dict.items():
                    f_val = [int(i.split('_')[-1]) for i in val]
                    ave_val = round(float(sum(f_val)) / float(len(f_val)), 2)
                    index1 = list(previous_strain_dict.keys()).index(key)
                    new_key = str(index1) + '_' + str(ave_val)
                    previous_strain_dict1[new_key] = val

                return previous_strain_dict1
            previous_strain_dict = previous_strain_dict
            key_list = sorted(previous_strain_dict.keys(), key=lambda x: len(previous_strain_dict[x]))
            key = key_list[ROUND + 1]
            val = previous_strain_dict[key]
            # print(previous_strain_dict.keys())

            ROUND += 1


# ######################### check R and R+1 ####################################################
#
# def check_R1_strain(n1_nodes_strain_dict, n1_nodes_vector_dict, BIC_file):
#     previous_strain_dict = deepcopy(n1_nodes_strain_dict)
#     f1 = open(BIC_file, 'w')
#
#     print(n1_nodes_strain_dict)
#
#     for key in sorted(n1_nodes_strain_dict.keys(), key=lambda x: int(x.split('_')[0])):
#         val = n1_nodes_strain_dict[key]
#         previous_input = {}
#         previous_input_nodes_vector_dict = {}
#         for i in val:
#             previous_input_nodes_vector_dict[i] = n1_nodes_vector_dict[i]
#         previous_input[key] = val
#         previous_result = get_likehood(previous_input, previous_input_nodes_vector_dict)
#         previous_likehood = previous_result[0]
#         previous_pi_list = previous_result[1]
#         previous_non_zero_list = previous_result[2]
#
#         k_mean_group = divide_into_two_group(previous_input_nodes_vector_dict)
#         new_group1 = k_mean_group[0]
#         new_group2 = k_mean_group[1]
#
#         # if len(new_group1) > 0.1 * len(val) and len(new_group2) > 0.1 * len(val):
#         if len(new_group1) > 5 and len(new_group2) > 5:
#             current_strain_dict = deepcopy(previous_strain_dict)
#
#             current_strain_dict[key + '_new1'] = new_group1
#             current_strain_dict[key + '_new2'] = new_group2
#
#             del current_strain_dict[key]
#
#             current_input = {}
#             current_input[key + '_new1'] = new_group1
#             current_input[key + '_new2'] = new_group2
#             current_input_nodes_vector_dict = previous_input_nodes_vector_dict
#             current_result = get_likehood(current_input, current_input_nodes_vector_dict)
#             current_likehood = current_result[0]
#             current_pi_list = current_result[1]
#
#             previous_parameter_num = get_parameter_number(previous_pi_list)
#
#             previous_nodes_num_list = [len(n1_nodes_strain_dict[key])]
#
#             previous_point_num = get_data_point_number1(previous_non_zero_list, previous_nodes_num_list)
#
#             previous_bic = BIC(previous_likehood, previous_parameter_num, previous_point_num)
#
#             current_parameter_num = get_parameter_number(current_pi_list)
#
#             current_point_num = previous_point_num
#             current_bic = BIC(current_likehood, current_parameter_num, current_point_num)
#
#             f1.write(str(previous_strain_dict.keys()) + '\t' + str(current_strain_dict.keys()) + '\t' + str(
#                 previous_bic) + '\t' + str(current_bic) + '\t' + str(
#                 [previous_likehood, previous_parameter_num, previous_point_num]) + '\t'
#                      + str([current_likehood, current_parameter_num, current_point_num])  + '\n')
#
#             diff = abs(previous_likehood - current_likehood)
#             deg_diff = abs(previous_parameter_num - current_parameter_num)
#
#             likelihood_p = get_likelihood_ratio(diff, deg_diff)
#
#             if likelihood_p < 1e-110:   ### 1e-110 is accepted
#                 previous_strain_dict = current_strain_dict
#                 # print(previous_strain_dict.keys())
#             else:
#                 previous_strain_dict = previous_strain_dict
#                 # print(previous_strain_dict.keys())
#         else:
#             previous_strain_dict = previous_strain_dict
#     f1.close()
#
#     previous_strain_dict1 = {}
#
#     for key, val in previous_strain_dict.items():
#         f_val = [int(i.split('_')[-1]) for i in val]
#         ave_val = round(float(sum(f_val)) / float(len(f_val)), 2)
#         new_key = key.split('_')[0] + '_' + str(ave_val)
#         previous_strain_dict1[new_key] = val
#
#     return previous_strain_dict1


def check_divide_into_two(key, nodes_list, n1_nodes_vector_dict):
    previous_input = {}
    previous_input_nodes_vector_dict = {}
    for i in sorted(nodes_list):
        previous_input_nodes_vector_dict[i] = n1_nodes_vector_dict[i]

    previous_input[key] = nodes_list
    previous_result = get_likehood(previous_input, previous_input_nodes_vector_dict)
    previous_likehood = previous_result[0]
    previous_pi_list = previous_result[1]
    print(key)
    k_mean_group = divide_into_two_group(previous_input_nodes_vector_dict)
    new_group1 = k_mean_group[0]
    new_group2 = k_mean_group[1]

    # if len(new_group1) > 0.1 * len(val) and len(new_group2) > 0.1 * len(val):
    if len(new_group1) > 5 and len(new_group2) > 5:
        current_input = {}
        current_input[key + '_new1'] = new_group1
        current_input[key + '_new2'] = new_group2
        current_input_nodes_vector_dict = previous_input_nodes_vector_dict
        current_result = get_likehood(current_input, current_input_nodes_vector_dict)
        current_likehood = current_result[0]
        current_pi_list = current_result[1]

        previous_parameter_num = get_parameter_number(previous_pi_list)

        current_parameter_num = get_parameter_number(current_pi_list)

        diff = abs(previous_likehood - current_likehood)
        deg_diff = abs(previous_parameter_num - current_parameter_num)

        likelihood_p = get_likelihood_ratio(diff, deg_diff)


        score1 = get_score1(diff,deg_diff)
        score2 = get_score2(diff,deg_diff)
        score3 = get_score3(diff,deg_diff)
        print('likelihood_p,def_diff')
        print(likelihood_p, deg_diff,score3)
        # if likelihood_p < 1e-170 or deg_diff==0:  ### 1e-110 is accepted

        if 0 < likelihood_p < 1e-110 or deg_diff == 0:  ### 1e-110 is accepted
            if score3 >= 3:
                return current_input
            else:
                return previous_input
            # return current_input
        elif likelihood_p == 0:
            if round(score3,1) >= 4.5:
                return current_input
            else:
                return previous_input
        else:
            return previous_input

    else:
        return previous_input


def check_divide_into_two1(divide_result, n1_nodes_vector_dict):
    current_divide_result = divide_result.copy()

    previous_key = list(set(current_divide_result.keys()))

    ROUND = 1
    while True:
        for key, val in current_divide_result.copy().items():
            nodes_list = val
            current_result = check_divide_into_two(key, nodes_list, n1_nodes_vector_dict)
            del current_divide_result[key]
            current_divide_result.update(current_result)

        current_key = list(set(current_divide_result))

        if sorted(previous_key) == sorted(current_key):
            return current_divide_result

        previous_key = current_key

        ROUND += 1

        print('devide_into_two1:ROUND-----------' + str(ROUND))
        print(current_divide_result.keys())


######################### check R and R+1 ####################################################

def check_R1_strain(n1_nodes_strain_dict, n1_nodes_vector_dict):
    # previous_strain_dict = deepcopy(n1_nodes_strain_dict)

    previous_strain_dict = {}

    for key in sorted(n1_nodes_strain_dict.keys(), key=lambda x: int(x.split('_')[0])):
        val = n1_nodes_strain_dict[key]
        previous_input = {}
        previous_input_nodes_vector_dict = {}
        for i in val:
            previous_input_nodes_vector_dict[i] = n1_nodes_vector_dict[i]
        previous_input[key] = val

        nodes_list = val
        divide_result = check_divide_into_two(key, nodes_list, n1_nodes_vector_dict)
        if len(divide_result) == 1:
            previous_strain_dict.update(divide_result)
        else:
            divide_result1 = check_divide_into_two1(divide_result, n1_nodes_vector_dict)
            previous_strain_dict.update(divide_result1)

    previous_strain_dict1 = {}

    print(previous_strain_dict.keys())

    for key, val in previous_strain_dict.items():
        f_val = [int(i.split('_')[-1]) for i in val]
        ave_val = round(float(sum(f_val)) / float(len(f_val)), 2)
        index1 = list(previous_strain_dict.keys()).index(key)
        new_key = str(index1) + '_' + str(ave_val)
        previous_strain_dict1[new_key] = val

    return previous_strain_dict1


# def compare_two_vector(vector1, vector2, cutoff):
#     label = []
#     non_zero = 0
#     for i in range(len(vector1)):
#         # if vector1[i] !=0 or vector2[i] != 0:
#             non_zero +=1
#             if vector1[i] <= vector2[i] + cutoff :
#                 label.append(1)
#
#     if len(label) >= non_zero*0.9:
#         return True
#     else:
#         return False
#
#
# def compare_two_vector_absolute(vector1, vector2, cutoff):
#     label = []
#     non_zero = 0
#     for i in range(len(vector1)):
#         # if vector1[i] != 0 or vector2[i] != 0:
#             non_zero +=1
#             if abs(vector1[i] - vector2[i]) <= cutoff:
#                 label.append(1)
#
#     # if len(label) >= len(vector1) - 2:
#     if len(label) >= non_zero*0.8:
#         return True
#     else:
#         return False


def compare_two_vector(vector1, vector2, cutoff):
    label = []
    non_zero = 0
    for i in range(len(vector1)):
        # if vector1[i] !=0 or vector2[i] != 0:
            non_zero +=1
            if vector1[i] <= vector2[i] + cutoff:
                label.append(1)
    print(vector1,vector2)
    print(len(label),non_zero,round(non_zero*0.9))
    if len(label) >= round(non_zero*0.9):
        return True
    else:
        return False


def compare_two_vector_absolute(vector1, vector2, cutoff):
    label = []
    non_zero = 0
    for i in range(len(vector1)):
        # if vector1[i] != 0 or vector2[i] != 0:
            non_zero += 1
            if round(abs(vector1[i] - vector2[i])) <= cutoff:
                label.append(1)
    # print(vector1,vector2)
    # if len(label) >= len(vector1) - 2:
    print(vector1,vector2)
    print(len(label),non_zero,round(0.8 * non_zero),non_zero - 2)
    if len(label) >= round(0.9 * non_zero) or (len(label) >= round(0.9 * non_zero) - 2):
    # if len(label) == len(vector1):
        return True
    else:
        return False


from itertools import combinations
import networkx as nx


def get_shared_nodes(R_strain_dict, nodes_vector_dict):
    R_strain_vector = get_R_strain_vector(R_strain_dict, nodes_vector_dict)

    ################################# step11_A select similar strain for each strain ###################################
    selected_strain_dict = {}
    for key, val in R_strain_vector.items():
        temp_list = []
        for key1, val1 in R_strain_vector.items():
            # if key != key1 and len(val) < len(val1):
            node_num = len(R_strain_dict[key])
            node_num1 = len(R_strain_dict[key1])
            if key != key1 and (2*node_num + 5)  < node_num1: ##### in order to keep the strains with less nodes
                print(key,key1,len(R_strain_dict[key]),len(R_strain_dict[key1]))
                label = compare_two_vector(val1, val, 1)  ############## compare whether val1 <= val +1

                if label == True:
                    temp_list.append(key1)

        selected_strain_dict[key] = temp_list

    print('step11_A')
    print(selected_strain_dict)

    ################################# step11_B  select Q pair strain for the selected strains###########################

    print('step11_B')
    selected_strain_pair_dict = {}

    for key, val in selected_strain_dict.items():
        print(key)
        pair_list = list(combinations(val, 2))

        selected_pair = []

        for pair in pair_list:

            key1 = pair[0]
            key2 = pair[1]

            val1 = [sum(i) for i in zip(R_strain_vector[key1], R_strain_vector[key2])]
            val2 = R_strain_vector[key]
            label = compare_two_vector(val1, val2, 1)

            if label == True:
                selected_pair.append(pair)

        selected_strain_pair_dict[key] = selected_pair
    print('step11_B')
    print(selected_strain_pair_dict)
    ################################## step8_C ######################

    assign_to_dict = {}
    for key, val in selected_strain_pair_dict.items():
        val2 = R_strain_vector[key]

        if len(val) == 1:

            key1 = val[0][0]
            key2 = val[0][1]

            val1 = [sum(i) for i in zip(R_strain_vector[key1], R_strain_vector[key2])]

            label = compare_two_vector_absolute(val1, val2, 1.5)  # this could be set to be 2
            diff = [round(abs(val2[i] - val1[i]), 3) for i in range(len(val1))]

            print('diff')
            print(diff)

            if label == True:
                assign_to_dict[key] = list(val[0])

        # if len(val) == 2:
        #     strain_list = list(set(list(val[0]) + list(val[1])))
        #
        #     vector_list = [R_strain_vector[i] for i in strain_list]
        #
        #     val1 = []
        #
        #     for i in range(len(vector_list[0])):
        #         temp = 0
        #         for j in range(len(vector_list)):
        #             temp += vector_list[j][i]
        #
        #         val1.append(temp)
        #     # print([round(i, 3) for i in val1])
        #     # print([round(i, 3) for i in val2])
        #     label = compare_two_vector_absolute(val1, val2, 1.5)  # this could set to be 2
        #
        #     diff = [round(abs(val2[i] - val1[i]), 3) for i in range(len(val1))]
        #
        #     if label == True:
        #         assign_to_dict[key] = list(set(list(val[0]) + list(val[1])))

        if len(val) > 1:
            G = nx.Graph()
            G.add_edges_from([i for i in val])
            clq_list = sorted(list(nx.clique.find_cliques(G)),key=lambda x: len(x),reverse= True)
            # clq_list = list(nx.clique.find_cliques(G))
            print('clq_list')
            print(clq_list)

            selected_clque_list = []
            diff_clique_list = []
            for clq in clq_list:
                strain_list = clq

                vector_list = [R_strain_vector[i] for i in strain_list]

                val1 = []

                for i in range(len(vector_list[0])):
                    temp = 0
                    for j in range(len(vector_list)):
                        temp += vector_list[j][i]

                    val1.append(temp)
                # print([round(i, 3) for i in val1])
                # print([round(i, 3) for i in val2])
                label = compare_two_vector_absolute(val1, val2, 1.5)  # this could set to be 2

                diff = [round(abs(val2[i] - val1[i]), 3) for i in range(len(val1))]

                if label == True:
                    print('test')
                    print(clq)
                    selected_clque_list.append(clq)
                    diff_clique_list.append(sum(diff))


            if len(diff_clique_list) > 1:
                if len(selected_clque_list[0]) == len(selected_clque_list[1]):
                    min_value = min(diff_clique_list)
                    min_index = diff_clique_list.index(min_value)
                    selected_clf = selected_clque_list[min_index]
                    assign_to_dict[key] = selected_clf

                else:
                    num_list = [len(i) for i in selected_clque_list]
                    max_value = max(num_list)
                    max_index = num_list.index(max_value)
                    assign_to_dict[key] = selected_clque_list[max_index]

            if len(diff_clique_list) == 1:
                assign_to_dict[key] = selected_clque_list[0]

    print('ste11_C')
    print(assign_to_dict)

    assign_to_dict1 = {}

    # for key,vals in assign_to_dict.items():
    #
    #     new_val = []
    #     for val in vals:
    #         if val not in list(assign_to_dict.keys()):
    #             new_val.append(val)
    #
    #     assign_to_dict1[key] = new_val
    if len(assign_to_dict.keys()) > 0:
        assigned_key = [list(assign_to_dict.keys())[0]]
        assigned_to_key = assign_to_dict[assigned_key[0]]

        for key, vals in assign_to_dict.items():
            if key in assigned_key:

                assign_to_dict1[key] = vals
            else:
                if key not in assigned_to_key:
                    new_value = []
                    for val in vals:
                        if val not in assigned_key:
                            new_value.append(val)

                    assign_to_dict1[key] = new_value

                    assigned_key.append(key)

    print(assign_to_dict1)
    for key, vals in assign_to_dict1.items():
        assign_nodes = R_strain_dict[key]
        for node in assign_nodes:
            for val in vals:
                R_strain_dict[val].append(node)
        del R_strain_dict[key]

    print(R_strain_dict.keys())

    return R_strain_dict


def check_strain_in_sample(strain_dict, nodes_vector_dict):
    Y_list = []

    sample_num = len(list(nodes_vector_dict.values())[0])
    strain_num = len(strain_dict.keys())
    for j in range(sample_num):
        y_temp = []
        for r in range(strain_num):
            r_nodes_list = set(list(strain_dict.values())[r])
            non_zero_list = []
            for r_nodes_list_temp in r_nodes_list:
                if nodes_vector_dict[r_nodes_list_temp][j] > 0:
                    non_zero_list.append(1)

            if len(non_zero_list) >= 5:
                y_temp.append(1)

            else:
                y_temp.append(0)

        Y_list.append(y_temp)

    return Y_list


def get_strain(merged_file, mixture_strain_file, frequency_file_loc, result_dir):
    print('###### getting nodes #######')
    n1_nodes = get_nodes(merged_file)
    print('###### getting mixture R strain #######')
    R_strain_dict = get_R_strain(mixture_strain_file)

    print('###### getting nodes dictionary #######')
    R_nodes = []
    for i in list(R_strain_dict.values()):
        R_nodes += i
    f1 = open(result_dir + '/MixtureS_result.bed', 'w')
    for key, val in R_strain_dict.items():
        f1.write('>strain_' + key + '\t' + str(val) + '\n')
    f1.close()
    #
    # ################ get_all_node and filter ########################
    #
    reads_file_list = os.listdir(frequency_file_loc)
    # # # # # # #
    filter_reads_file = frequency_file_loc

    #
    # #############################################################
    #
    dict_all = get_nodes_frequency_vectors(n1_nodes, R_nodes, filter_reads_file)
    n1_nodes_vector_dict = dict_all[0]

    f1_1 = open(result_dir + '/n1_nodes_vector_normalized.bed', 'w')
    for key, val in n1_nodes_vector_dict.items():
        f1_1.write(key + '\t' + str(val) + '\n')
    f1_1.close()
    #
    n1_nodes_vector_dict_original = dict_all[2]
    #
    f1_2 = open(result_dir + '/n1_nodes_vector_original.bed', 'w')
    for key, val in n1_nodes_vector_dict_original.items():
        f1_2.write(key + '\t' + str(val) + '\n')
    f1_2.close()

    R_nodes_vector_dict = dict_all[1]
    # print(n1_nodes_vector_dict)
    # print(R_nodes_vector_dict)
    print('###### getting mixture R strain vector #######')
    R_strain_vector_dict = get_R_strain_vector(R_strain_dict, R_nodes_vector_dict)

    ################################### step8 ############################################

    print('step8')
    f1 = result_dir + '/round_nums_step8'
    n1_nodes_strain_dict = get_n1_nodes_strain(n1_nodes_vector_dict, R_strain_vector_dict, f1)

    f2 = open(result_dir + '/R_strain_step8.bed', 'w')
    for key, val in n1_nodes_strain_dict.items():
        # print(key, val)
        f2.write('>strain_' + key + '\t' + str(val) + '\n')
    f2.close()

    ############################################## step9 check R and R +1 #######################

    print('step9_divide_to_two_group')

    previous_strain_dict = check_R1_strain(n1_nodes_strain_dict, n1_nodes_vector_dict)

    print(previous_strain_dict.keys())

    ########################################### step9 assign nodes ##########################

    n1_strain_dict = deepcopy(previous_strain_dict)

    R_strain_vector_dict = get_R_strain_vector(n1_strain_dict, n1_nodes_vector_dict)
    # print(R_strain_vector_dict)
    print('step9_assign_nodes')
    f1 = result_dir + '/round_nums_step9'
    n1_nodes_strain_dict = get_n1_nodes_strain(n1_nodes_vector_dict, R_strain_vector_dict, f1)

    ############################################ step9 check R and R-1 strain #################
    BIC_9 = result_dir + '/R_strain_step9_test_BIC.bed'

    print('step9')

    previous_strain_dict = check_R_strain(n1_nodes_strain_dict, n1_nodes_vector_dict, BIC_9)

    f3 = open(result_dir + '/R_strain_step9.bed', 'w')
    for key, val in previous_strain_dict.items():
        f3.write('>strain_' + key + '\t' + str(val) + '\n')
    f3.close()

    ################ step10 assign the shared nodes ############################

    final_strain_dict = get_shared_nodes(previous_strain_dict, n1_nodes_vector_dict)

    f5 = open(result_dir + '/R_strain_step10.bed', 'w')
    for key, val in final_strain_dict.items():
        f5.write('>strain_' + key + '\t' + str(val) + '\n')

    f5.close()

    ############### assign the unnomal nodes again ###############################

    strain_dict_ocurs = {}
    for key, val in final_strain_dict.items():
        X1 = [n1_nodes_vector_dict[i] for i in val]

        X2 = []

        for i in X1:

            temp = []

            for item in i:
                if int(item) > 0:
                    temp.append(1)
                else:
                    temp.append(0)
            X2.append(temp)

        strain_dict_ocurs[key] = X2

    nodes_num = sum([len(i) for i in list(strain_dict_ocurs.values())])

    # nodes_num = len(n1_nodes)

    assinged_key = []

    cutoff_list = []

    for key, val in strain_dict_ocurs.items():

        sample_size = len(val[0])

        sum_list = []
        for i in range(sample_size):
            temp = []
            for j in val:
                temp.append(j[i])

            temp1 = sum(temp)
            sum_list.append(temp1)

        sum_list1 = []

        for i in sum_list:
            if i <= 5:
                sum_list1.append(0)
            else:
                sum_list1.append(i)

        ave = [i / nodes_num for i in sum_list1]

        # cutoff = sum(ave1) / len(ave1)

        cutoff = max(ave)

        # print('cutoff')
        # print(cutoff)

        cutoff_list.append(cutoff)

    cutoff_list_100 = normalized_to_max1(cutoff_list,100)
    # print(cutoff_list_100)

    for i in range(len(cutoff_list_100)):
            if cutoff_list_100[i] <= 12:
                assigned_key2 = list(strain_dict_ocurs.keys())[i]
                val_dict = {}
                val = final_strain_dict[assigned_key2]
                for i in val:
                    val_dict[i] = n1_nodes_vector_dict[i]
                current_strain_dict = deepcopy(final_strain_dict)
                del current_strain_dict[assigned_key2]

                strain_vector_dict = get_R_strain_vector(current_strain_dict, n1_nodes_vector_dict)
                final_strain_dict = get_nodes_assigned(val_dict, strain_vector_dict, current_strain_dict)

                print('assign unnormal nodes')
                print(final_strain_dict.keys())

    ################## get the original coverage of strain ##########################
    all_nodes_in_strain = []
    print(final_strain_dict.keys())
    for i in list(final_strain_dict.values()):
        all_nodes_in_strain += i

    nodes_number_dict = Counter(all_nodes_in_strain)
    f5_1 = open(result_dir + '/MSMS_result.bed', 'w')

    final_strain_dict_original = {}
    for key, val in final_strain_dict.items():
        val1 = [i for i in val if nodes_number_dict[i] == 1]

        fre_num = [sum(n1_nodes_vector_dict_original[i]) for i in val1]

        ave_fre = round(float(sum(fre_num)) / float(len(fre_num)), 2)

        index1 = list(final_strain_dict.keys()).index(key)

        final_strain_dict_original['strain_' + str(index1) + '_' + str(ave_fre)] = val

        f5_1.write('>strain_' + str(ave_fre) + '\t' + str(val) + '\n')

    f5_1.close()

    ############### get_strain_in_each_sample #############################

    strain_in_sample_list = check_strain_in_sample(final_strain_dict_original, n1_nodes_vector_dict_original)

    sample_name = reads_file_list
    strain_name = list(final_strain_dict_original.keys())

    f6 = open(result_dir + '/strain_in_each_sample.txt', 'w')

    f6.write('sample_name' + '\t' + '\t'.join(strain_name) + '\n')

    # print(sample_name)

    for i in range(len(strain_in_sample_list)):
        f6.write(sample_name[i] + '\t' + '\t'.join([str(item) for item in strain_in_sample_list[i]]) + '\n')

    f6.close()




















