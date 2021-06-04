from os.path import isfile, join
from os import listdir
from rpy2.robjects.packages import importr

rStats = importr('stats')
import rpy2.robjects as robjects
from operator import itemgetter
import numpy as np
from copy import deepcopy
import os.path
import jenkspy
from scipy.stats import chi2
from multiprocessing import Pool
from collections import defaultdict
import sys

CRITERIA = 1e-6
polymorphic_cutoff = 3
REF_COUNT_NORMALIZED = 5
SECOND_COUNT_NORMALIZED = 5
PI_THRESHOLD = 0.05
ALPHA_THRESHOLD = 0.05
################################################################################
'''
multinomial distribution:

(1)R:
> x<-c(1,2,3)
> prob<-c(0.2,0.3,0.5)
> dmultinom(x,prob=prob)
[1] 0.135
-------------------------------
> x<-c(20,40,40)
> dmultinom(x,prob=prob)
[1] 0.0006680965


(2)Python:
multinomial.pmf([20,40,40],n=100,p=[0.2,0.3,0.5])
Out[35]: 0.000668096469909772
multinomial.pmf([1,2,3],n=6,p=[0.2,0.3,0.5])
Out[36]: 0.13500000000000018

In this script, we don't consider the ACGT position. We only consider the 
proportion of each 

(1) x2 / (x1 + x2)
(2) binomial distribution
(3) n -> n + 1: stop at if now_numHaps == previous_numHaps or now_numHaps < current_predict




EM_V6: implement Likelihood-ratio test
EM_v7: implement BIC
EM_v11: updated the one haplotype:
(1) 
'''


################################################################################
def normalized_to_100(count):
    freq = []
    coe = 100.0 / sum(count)
    for i in range(len(count)):
        freq.append(int(round(count[i] * coe, 0)))

    # make it is A+C+G+T=100, may be not necessay
    if sum(freq) != 100:
        difference = 100 - sum(freq)
        maxIndex = [index for index, value in enumerate(freq) if
                    value == max(freq)]
        freq[maxIndex[0]] += difference

    return freq


def normalized_to_100_frequency(count):
    freq = []
    coe = 100.0 / sum(count)
    for i in range(len(count)):
        freq.append(int(round(count[i] * coe, 0)))

    # make it is A+C+G+T=100, may be not necessay
    if sum(freq) != 100:
        difference = 100 - sum(freq)
        maxIndex = [index for index, value in enumerate(freq) if
                    value == max(freq)]
        freq[maxIndex[0]] += difference

    freq = [round(freq[i] / 100.0, 2) for i in range(len(freq))]

    return freq


# we do not consider the sequence error here, so initial the parameter based on
# quantile for all haplotypes.
def initial_parameter_position_no_error_v2(polyCounts, numHaps, data_ref):
    polyCounts.sort(key=lambda x: x[2])
    # method3: Jenks natural breaks, optimization
    # order_2nd = [item[2] for item in polyCounts]
    order = [item[2] * 1.0 / (item[2] + item[1]) for item in polyCounts]
    breaks = jenkspy.jenks_breaks(order, nb_class=numHaps)
    breaks.sort()
    cluster = [[] for _ in range(numHaps)]
    count = [0] * numHaps

    for item in order:

        for i in range(1, len(breaks)):
            if i == len(breaks) - 1:
                if item >= breaks[i - 1] and item <= breaks[i]:
                    cluster[i - 1].append(item)
                    count[i - 1] += 1
                    break
            else:
                if item >= breaks[i - 1] and item < breaks[i]:
                    cluster[i - 1].append(item)
                    count[i - 1] += 1
                    break
    median_list = []
    for item in cluster:
        if len(item) == 0:
            median_list.append(0.0)
        else:
            median_list.append(np.median(item))

    # just mark sure sum of pi list = 1
    pi = [item * 1.0 / sum(median_list) for item in median_list]
    pi[-1] = 1.0 - sum(pi[:-1])

    alpha = []
    for i in range(len(count)):
        alpha.append(count[i] * 1.0 / sum(count))
    alpha[-1] = 1.0 - sum(alpha[:-1])
    return pi, alpha


def read_polymorphicsSite_label(inName):
    data = []
    with open('%s' % inName) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            line = [int(one) for one in line[:-1]] + [line[-1]]
            data.append(line)
    return data


def read_polymorphicsSite(inName):
    data = []
    with open('%s' % inName) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            line = [int(one) for one in line]
            data.append(line)
    return data


def compare_with_trueLabel(polyCounts_cluster, data, resDir):
    # get all labels from
    predict_labels = [value for key, value in polyCounts_cluster.items()]
    true_labels = []
    for item in data:
        tmp = []
        if ';' in item[5]:
            tmp = [int(one) for one in item[5].split(';')]
        else:
            tmp.append(int(item[5]))
    # true_labels = [item[5] for item in data]
    labels_all_true = []
    for item in true_labels:
        labels_all_true += item

    labels_all = list(set(labels_all_true + predict_labels))

    contingency_individual = {}
    contingency_total = {'TRUE': 0, 'FALSE': 0, 'TOTAL': 0}
    for label in labels_all:
        contingency_individual.update(
            {label: {'TRUE': 0, 'FALSE': 0, 'TOTAL': 0}})

    for i in range(len(data)):
        loc = data[i][0]
        true_label = []
        if ';' in data[i][5]:
            true_label = [int(one) for one in data[i][5].split(';')]
        else:
            true_label = [int(data[i][5])]
        # true_label = data[i][5]
        predict_label = polyCounts_cluster[loc]

        # if true_label==predict_label:
        if predict_label in true_label:
            contingency_individual[predict_label]['TRUE'] += 1
            contingency_total['TRUE'] += 1
        else:
            contingency_individual[predict_label]['FALSE'] += 1
            contingency_total['FALSE'] += 1

        contingency_individual[predict_label]['TOTAL'] += 1
        contingency_total['TOTAL'] += 1

    # write into files
    order = ['TRUE', 'FALSE', 'TOTAL']
    resName = resDir + '/polymorphicSites_statistic'
    with open('%s' % resName, 'w') as f:
        for label in labels_all:
            tmp = '\t'.join([str(label)] +
                            [str(contingency_individual[label][key]) for key in
                             order])
            f.write('%s\n' % tmp)

        tmp = '\t'.join(
            ['total'] + [str(contingency_total[key]) for key in order])
        f.write('%s' % tmp)

    # write true and predict label into file
    resName = resDir + '/polymorphicSites'
    with open('%s' % resName, 'w') as f:
        for i in range(len(data)):
            loc = data[i][0]
            true_label = data[i][5]
            predict_label = polyCounts_cluster[loc]
            if true_label == predict_label:
                label = 'True'
            else:
                label = 'False'

            tmp = data[i] + [label, predict_label]
            tmp = '\t'.join([str(one) for one in tmp])
            f.write('%s\n' % tmp)


def prediction(data, numHaps, proportion, alpha):
    polyCounts = deepcopy(data)
    labels_all = [0] + [i + 1 for i in range(numHaps)]

    # print parameter
    print('\t'.join([''] + [str(one) for one in labels_all]))
    print('\t'.join(['pi'] + [str(one) for one in proportion]))
    for index, value in enumerate(labels_all):
        print('\t'.join([str(value)] + [str(one) for one in alpha[index]]))

    # normalized polymorphic sites to sum of 100. A+C+G+T=100
    for i in range(len(polyCounts)):
        polyCounts[i][1:5] = normalized_to_100(polyCounts[i][1:5])
        # print 'here'

    # two way sort
    for i in range(len(polyCounts)):
        polyCounts[i][1:5] = sorted(polyCounts[i][1:5], reverse=True)
    # polyCounts.sort(key=itemgetter(1),reverse=True)
    polyCounts.sort(key=lambda item: (item[1], item[2], item[3], item[4]),
                    reverse=True)

    # calculate the posterior
    # E-Step: calculate the posterior given observed data
    # (1)likelihood matrix for observed data given haplotypes, p(x|z)
    likelihood = []
    for dataID in range(len(polyCounts)):
        item = polyCounts[dataID]
        probList = []
        for hapID in range(numHaps + 1):
            alphaList = alpha[hapID]
            dataList = polyCounts[dataID][1:3]

            # #python
            # prob = multinomial.pmf(item[1:], n=sum(item[1:]), p=likelihood_parameter[hapID])
            # R
            x = robjects.FloatVector(dataList)
            p = robjects.FloatVector(alphaList)
            prob = rStats.dmultinom(x=x, prob=p)[0]
            probList.append(prob)
        likelihood.append(probList)

    # (2)calculate the posterior given data x, p(z|x). This is for polymorphic sites
    posterior = []
    for dataID in range(len(polyCounts)):
        probList = []

        denominator = 0.0
        for hapID in range(numHaps + 1):
            prob = likelihood[dataID][hapID] * proportion[hapID]
            denominator += prob

        for hapID in range(numHaps + 1):
            prob = likelihood[dataID][hapID] * proportion[hapID] / denominator
            probList.append(prob)

        posterior.append(probList)

    # prediction the label according the posterior
    labels = []
    for i in range(len(posterior)):
        maxIndex = \
        max((value, index) for index, value in enumerate(posterior[i]))[1]
        label = maxIndex
        labels.append(label)

    return polyCounts, labels


# similar to normalized_filter(), but not filter anything. just format the data.
def format_data(data, data_ref):
    data_raw = deepcopy(data)

    # normalized polymorphic sites to sum of 100. A+C+G+T=100. and only pick
    # the data in the format = [ref, largest_execept_ref]
    polyCounts = []
    orderMap = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    count_third = 0
    for i in range(len(data_raw)):
        # tmp = normalized_to_100(data_raw[i][1:5])
        tmp = data_raw[i][1:5]

        # get ref
        loc = data_raw[i][0]
        ref_index = orderMap[data_ref[loc]]

        ref_count = tmp[ref_index]
        # #filter ref count
        # if ref_count < REF_COUNT_NORMALIZED:
        #     continue

        # the largest number in each row except reference.
        second_count = 0
        second_index = 0
        for index in range(4):
            if index == ref_index:
                continue
            if tmp[index] > second_count:
                second_count = tmp[index]
                second_index = index

        # third letter larger than ref
        third_count = 0
        for index in range(4):
            if index == second_index or index == ref_index:
                continue
            if tmp[index] > ref_count:
                third_count = tmp[index]

        if third_count == 0:
            one = ref_count
            two = second_count
        else:
            one = second_count
            two = third_count
            count_third += 1

        # if one <= REF_COUNT_NORMALIZED or two <= SECOND_COUNT_NORMALIZED:
        #     continue
        tmp = [data_raw[i][0], one, two] + tmp

        polyCounts.append(tmp[:])

    # polyCounts.sort(key=itemgetter(1),reverse=True)
    polyCounts.sort(key=lambda item: (item[1], item[2]), reverse=True)

    return polyCounts


def normalized_filter(data, data_ref):
    data_raw = deepcopy(data)
    # print(data_raw)
    # normalized polymorphic sites to sum of 100. A+C+G+T=100. and only pick
    # the data in the format = [ref, largest_execept_ref]
    polyCounts = []
    orderMap = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    count_third = 0
    for i in range(len(data_raw)):
        tmp = normalized_to_100(data_raw[i][1:5])

        # get ref
        loc = data_raw[i][0]
        ref_index = orderMap[data_ref[loc]]

        ref_count = tmp[ref_index]
        # filter ref count
        if ref_count < REF_COUNT_NORMALIZED:
            continue

        # the largest number in each row except reference.
        second_count = 0
        second_index = 0
        for index in range(4):
            if index == ref_index:
                continue
            if tmp[index] > second_count:
                second_count = tmp[index]
                second_index = index

        # third letter larger than ref
        third_count = 0
        for index in range(4):
            if index == second_index or index == ref_index:
                continue
            if tmp[index] > ref_count:
                third_count = tmp[index]

        if third_count == 0:
            one = ref_count
            two = second_count
        else:
            one = second_count
            two = third_count
            count_third += 1

        if one <= REF_COUNT_NORMALIZED or two <= SECOND_COUNT_NORMALIZED:
            continue
        tmp = [data_raw[i][0], one, two] + tmp

        polyCounts.append(tmp[:])

    print(count_third)
    # print(data_raw)
    # print(polyCounts)
    tmp = [len(data_raw), len(polyCounts),
           round(len(polyCounts) * 1.0 / len(data_raw), 5)]
    print('survived polymorphic sites after fitler each letter less than ' + \
          str(REF_COUNT_NORMALIZED) + ',' + str(
        SECOND_COUNT_NORMALIZED) + ': ' + \
          ' '.join([str(item) for item in tmp]))

    # polyCounts.sort(key=itemgetter(1),reverse=True)
    polyCounts.sort(key=lambda item: (item[1], item[2]), reverse=True)
    return polyCounts


# return l1 and l2 are independent?
def likelihood_ratio_test(likelihood_log, index):
    if index <= 2:
        return True

    l2 = likelihood_log[index]
    l1 = likelihood_log[index - 1]

    # LR = -2.0 * np.log(l1 / l2)
    LR = 2 * (l2 - l1)
    p = chi2.sf(LR, 1)
    '''
    null hypothesis: l1 and l2 are independent
    large LR -> small p-value from chi-squre(1) -> not reject null hypothesis
    -> l1 and l2 are independent
    small LR -> large p-value from chi-squre(1) -> reject null hypothesis
    -> l1 and l2 are not independent
    '''

    if p < 0.05:
        return True
    else:
        return False


def AIC(likelihood_log, index):
    if index <= 2:
        return True

    l2 = likelihood_log[index]
    l1 = likelihood_log[index - 1]

    # AIC = 2k - 2 logL
    aic1 = 2 * (2 * (index - 1)) - 2 * l1
    aic2 = 2 * (2 * index) - 2 * l2

    if aic1 <= aic2:
        return False
    else:
        return True


def BIC(likelihood_log, index, nums_data):
    if index <= 2:
        l2 = likelihood_log[index]
        bic2 = np.log(nums_data) * (2 * (index)) - 2 * l2
        return True, bic2

    l2 = likelihood_log[index]
    l1 = likelihood_log[index - 1]

    bic1 = np.log(nums_data) * (2 * (index - 1)) - 2 * l1
    bic2 = np.log(nums_data) * (2 * (index)) - 2 * l2

    if bic2 < bic1:
        return True, bic2
    else:
        return False, bic2


def calculate_BIC_one_hap(inFileName, genomeFileLoc):
    # read polymorphic sites and filter
    data = read_polymorphicsSite(inFileName)
    count_before = len(data)
    data = filter(data)

    count_after = len(data)
    tmp = [count_before, count_after,
           round(count_after * 1.0 / count_before, 5)]
    print(inFileName + ' ratio of polymorphic sites survived after filter 10% of total reads: ' + \
          ' '.join([str(item) for item in tmp]))

    # label polymorphic sites what is reference in that pos.
    polyCounts_refLabel = check_refernce_genome(data, genomeFileLoc)

    # normalized and filter individual letter
    data_normalized = normalized_filter(data, polyCounts_refLabel)

    order = [item[2] * 1.0 / (item[2] + item[1]) for item in data_normalized]
    pi = np.median(order)
    alpha = 1

    # log likelihood sum
    '''
    L(theta) = f(x1;theta)*f(x2;theta)*...*f(xN;theta)
    f(x1;theta) = lambda1*f(x1;theta1) + lambda2*f(x1;theta2) + ...+ lambdaK*f(x1;thetaK)
    N: samples, K: different distribution.

    logL(theta) = sum(logf(xi;theta))
    '''
    data = deepcopy(data_normalized)
    likelihood_sum = 0.00
    for dataID in range(len(data)):
        item = data[dataID]
        data_alpha_list = [item[1], item[2]]
        pi_list = [1.0 - pi, pi]

        x = robjects.FloatVector(data_alpha_list)
        p = robjects.FloatVector(pi_list)
        prob = rStats.dmultinom(x=x, prob=p)[0]

        # f(x1;theta)
        denominator = prob * alpha
        likelihood_sum += np.log(denominator)

    bic = np.log(len(data_normalized)) * (2 * (1)) - 2 * likelihood_sum

    return bic


def filter_one_hap(inFileName, genomeFileLoc):
    # read polymorphic sites and filter
    data = read_polymorphicsSite(inFileName)
    count_before = len(data)
    data = filter(data)

    count_after = len(data)
    tmp = [count_before, count_after,
           round(count_after * 1.0 / count_before, 5)]
    print(inFileName + ' ratio of polymorphic sites survived after filter 10% of total reads: ' + \
          ' '.join([str(item) for item in tmp]))

    # label polymorphic sites what is reference in that pos.
    polyCounts_refLabel = check_refernce_genome(data, genomeFileLoc)

    # normalized and filter individual letter
    data_normalized = normalized_filter(data, polyCounts_refLabel)

    # median_value = np.median([item[2] for item in data_normalized])
    snp_output = generate_final_SNPs_single_hap(data_normalized, polyCounts_refLabel)
    # if len(data_normalized) < 50 or median_value <= 10:
    if len(data_normalized) < 50:
        return True, snp_output

    return False, snp_output


def EM(data, numHaps, resDir, data_ref):
    # pi, likelihood_parameter = initial_parameter_position_new(polyCounts, numHaps,
    #                                                       data_ref)
    pi, likelihood_parameter = initial_parameter_position_no_error_v2(data,
                                                                      numHaps,
                                                                      data_ref)

    # write initialize info file
    resName = resDir + '/EM_initialize_pi'
    with open('%s' % resName, 'w') as f:
        tmp = '\t'.join([str(one) for one in pi])
        f.write('%s' % tmp)
    resName = resDir + '/EM_initialize_alpha'
    with open('%s' % resName, 'w') as f:
        for item in likelihood_parameter:
            # tmp = '\t'.join([str(one) for one in item])
            tmp = str(item)
            f.write('%s\n' % tmp)

    print('\t'.join([str(one) for one in pi]))
    for item in likelihood_parameter:
        # print '\t'.join([str(one) for one in item])
        print(str(item))

    resName = resDir + '/input_normalized'
    with open('%s' % resName, 'w') as f:
        for line in data:
            tmp = '\t'.join([str(one) for one in line])
            f.write('%s\n' % tmp)

    # (1)Changed to double sort format
    previous_pi = ''
    previous_likelihood = ''

    # write EM intermedia result into file
    resName = resDir + '/EM_intermedia_parameter'
    f_intermedia = open('%s' % resName, 'w')
    # tmp = '\t'.join([str(one) for one in pi])
    # f.write('%s' %tmp)

    ROUND = 1
    # EM algorithm
    while True:
        # E-Step: calculate the posterior given observed data
        # (1)likelihood matrix for observed data given haplotypes, p(x|z)
        likelihood = []

        # likelihood for polymorphic sites
        for dataID in range(len(data)):
            item = data[dataID]
            probList = []
            for hapID in range(numHaps):
                data_alpha_list = [item[1], item[2]]
                # alphaList = likelihood_parameter[hapID]
                pi_list = [1.0 - pi[hapID], pi[hapID]]
                # data_alpha_list = [100 - data_alpha, data_alpha]
                # data_alpha_list = [item[1], item[2]]

                # #python
                # prob = multinomial.pmf(item[1:], n=sum(item[1:]),
                # p=likelihood_parameter[hapID])
                # R
                x = robjects.FloatVector(data_alpha_list)
                p = robjects.FloatVector(pi_list)
                prob = rStats.dmultinom(x=x, prob=p)[0]
                probList.append(prob)
            likelihood.append(probList)

        # (2)calculate the posterior given data x, p(z|x). This is for
        # polymorphic sites
        posterior = []
        for dataID in range(len(data)):
            probList = []

            denominator = 0.0
            for hapID in range(numHaps):
                prob = likelihood[dataID][hapID] * likelihood_parameter[hapID]
                denominator += prob

            # if denominator == 0.0:
            #     denominator = 1e-10
            for hapID in range(numHaps):
                prob = likelihood[dataID][hapID] * likelihood_parameter[
                    hapID] / denominator
                probList.append(prob)

            posterior.append(probList)

        # M-step
        # (1)update alpha
        for hapID in range(numHaps):
            # pi[hapID] = sum([item[hapID] for item in posterior])/len(posterior)
            likelihood_parameter[hapID] = sum(
                [item[hapID] for item in posterior]) / len(
                posterior)

        # update pi
        for hapID in range(numHaps):

            W = []
            for dataID in range(len(data)):
                item = data[dataID]
                W.append(posterior[dataID][hapID] * 1.0 * item[2] / (
                            item[1] + item[2]))
            pi[hapID] = sum(W) / sum([item[hapID] for item in posterior])

        # check when the EM is stopped
        flag_converge = True
        if ROUND == 1:
            flag_converge = False
            previous_likelihood = deepcopy(likelihood_parameter)
            previous_pi = deepcopy(pi)
        else:
            # check pi
            for i in range(len(pi)):
                if abs(previous_pi[i] - pi[i]) > CRITERIA:
                    flag_converge = False

            previous_likelihood = deepcopy(likelihood_parameter)
            previous_pi = deepcopy(pi)

        if ROUND > 500:
            resTmp = resDir + '_converge_log'
            with open('%s' % resTmp, 'w') as f:
                f.write('True')
            flag_converge = True

        # print '#########################'
        print('EM round\t' + str(ROUND))
        print('\t'.join(['pi'] + [str(one) for one in pi]))
        print('\t'.join(['alpha'] + [str(one) for one in likelihood_parameter]))

        # write intermediate result into file
        f_intermedia.write('EM round\t' + str(ROUND))
        f_intermedia.write('\n')
        f_intermedia.write('\t'.join(['pi'] + [str(one) for one in pi]))
        f_intermedia.write('\n')
        f_intermedia.write(
            '\t'.join(['alpha'] + [str(one) for one in likelihood_parameter]))
        f_intermedia.write('\n')

        if flag_converge:
            f_intermedia.close()
            return ROUND, pi, likelihood_parameter

        ROUND += 1


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


def posterior_cluster(data_raw, pi, likelihood_parameter, numHaps, data_ref):
    data = deepcopy(data_raw)

    likelihood = []

    # likelihood for polymorphic sites
    for dataID in range(len(data)):
        item = data[dataID]
        probList = []
        for hapID in range(numHaps):
            data_alpha_list = [item[1], item[2]]
            # alphaList = likelihood_parameter[hapID]
            pi_list = [1.0 - pi[hapID], pi[hapID]]
            # data_alpha_list = [100 - data_alpha, data_alpha]
            # data_alpha_list = [item[1], item[2]]

            # #python
            # prob = multinomial.pmf(item[1:], n=sum(item[1:]),
            # p=likelihood_parameter[hapID])
            # R
            x = robjects.FloatVector(data_alpha_list)
            p = robjects.FloatVector(pi_list)
            prob = rStats.dmultinom(x=x, prob=p)[0]
            probList.append(prob)
        likelihood.append(probList)

    # log likelihood sum
    '''
    L(theta) = f(x1;theta)*f(x2;theta)*...*f(xN;theta)
    f(x1;theta) = lambda1*f(x1;theta1) + lambda2*f(x1;theta2) + ...+ lambdaK*f(x1;thetaK)
    N: samples, K: different distribution.

    logL(theta) = sum(logf(xi;theta))
    '''
    likelihood_sum = 0.00
    for dataID in range(len(data)):

        # f(x1;theta)
        denominator = 0.0
        for hapID in range(numHaps):
            prob = likelihood[dataID][hapID] * likelihood_parameter[hapID]
            denominator += prob
        likelihood_sum += np.log(denominator)

    # likelihood_sum = 0.00
    # for dataID in range(len(data)):
    #
    #     #f(x1;theta)
    #     denominator = 0.0
    #     for hapID in range(numHaps):
    #         prob = likelihood[dataID][hapID] * likelihood_parameter[hapID]
    #         denominator += prob
    #     likelihood_sum += np.log(denominator)

    # (2)calculate the posterior given data x, p(z|x). This is for
    # polymorphic sites
    posterior = []
    for dataID in range(len(data)):
        probList = []

        denominator = 0.0
        for hapID in range(numHaps):
            prob = likelihood[dataID][hapID] * likelihood_parameter[hapID]
            denominator += prob

        # if denominator == 0.0:
        #     denominator = 1e-10
        for hapID in range(numHaps):
            prob = likelihood[dataID][hapID] * likelihood_parameter[
                hapID] / denominator
            probList.append(prob)

        posterior.append(probList)
    print('posterior')
    print(posterior)
    # prediction the label according the posterior for sure label, not_sure=-1
    labels = []
    for i in range(len(posterior)):
        label = -1
        for index in range(len(posterior[i])):
            if posterior[i][index] >= 0.6:
                label = index
                break
        labels.append(label)
    print('labels')
    print(labels)
    #assign position to corresponding cluster
    cluster = defaultdict(list)
    for i in range(len(data)):
        if labels[i] != -1:
            cluster[data[i][0]].append(labels[i])
            continue

        # for index in range(len(posterior[i])):
        #     if posterior[i][index] >= 0.05:
        #         cluster[data[i][0]].append(index)

        max_value = max(posterior[i])
        index = posterior[i].index(max_value)
        cluster[data[i][0]].append(index)
    # cluster = defaultdict(list)
    #
    # for i in range(len(data)):
    #     max_value = max(posterior[i])
    #     index = posterior[i].index(max_value)
    #     cluster[data[i][0]].append(index)

    print('cluster:')
    print(numHaps)
    print(cluster)

    return cluster, likelihood_sum


def calculate_proportion(polyCounts, polyCounts_normalized, refLable,
                         clusterLabel, resDir):
    # step1: cluster the polymorphic sites into groups
    totalLabel = set()
    for key, value in clusterLabel.items():
        for one in value:
            totalLabel.add(one)
    totalLabel = list(totalLabel)

    # initial
    poly_cluster = {}
    poly_cluster_normalized = {}
    for label in totalLabel:
        poly_cluster.update({label: []})
        poly_cluster_normalized.update({label: []})

    orderMap = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    # cluster by normalized frequency
    for item in polyCounts_normalized:
        pos = item[0]
        # cov = item[2]
        cluster = clusterLabel[pos]
        for hapID in cluster:
            poly_cluster_normalized[hapID].append(
                1.0 * item[2] / (item[1] + item[2]))

    # cluster by the abosolute count
    for item in polyCounts:
        pos = item[0]
        # cov = item[1]
        if pos not in clusterLabel:
            continue
        cluster = clusterLabel[pos]
        for hapID in cluster:
            poly_cluster[hapID].append(1.0 * item[2] / (item[1] + item[2]))

    # not normalized to 100
    print('################')
    print('Proportions after normalized:')
    proportions_normalized = []
    for key, value in poly_cluster_normalized.items():
        print('\t'.join([str(key), str(np.median(value))]))
        proportions_normalized.append([str(key), np.median(value)])
    print('#######################')
    print('Proportions for the raw read counts:')
    proportions = []
    for key, value in poly_cluster.items():
        print('\t'.join([str(key), str(np.median(value))]))
        proportions.append([str(key), np.median(value)])
    print('#######################')

    # normalized to 100
    sum_raw = sum([item[1] for item in proportions])
    sum_normalized = sum([item[1] for item in proportions_normalized])

    for item in proportions:
        tmp = item[1] * 1.0 / sum_raw
        item.append(tmp)
    for item in proportions_normalized:
        tmp = item[1] * 1.0 / sum_normalized
        item.append(tmp)

    # normalized to 100
    for item in proportions:
        tmp = '\t'.join(['raw_reads:'] + [str(one) for one in item])
        print(tmp)
    for item in proportions_normalized:
        tmp = '\t'.join(['normalized:'] + [str(one) for one in item])
        print(tmp)

    resName = resDir + '/proportion'
    with open('%s' % resName, 'w') as f:
        for item in proportions:
            tmp = '\t'.join(['raw_reads:'] + [str(one) for one in item])
            f.write('%s\n' % tmp)

        for item in proportions_normalized:
            tmp = '\t'.join(['normalized:'] + [str(one) for one in item])
            f.write('%s\n' % tmp)

    return proportions_normalized, proportions


def decide_num_haplotype(pi, likelihood_parameter):
    count_similar = len(pi)
    similar_index = []
    for i in range(1, len(pi)):
        if pi[i] - pi[i - 1] <= PI_THRESHOLD:
            count_similar -= 1
            similar_index.append(i)
            similar_index.append(i - 1)

    for i in range(len(likelihood_parameter)):
        if likelihood_parameter[i] <= ALPHA_THRESHOLD:
            if i not in similar_index:
                count_similar -= 1

    return count_similar

# def decide_num_haplotype(pi, likelihood_parameter):
#     count_similar = len(pi)
#     similar_index = []
#     for i in range(len(pi)-1):
#         for j in range(i+1,len(pi)):
#             if pi[j] - pi[i] <= PI_THRESHOLD:
#                 count_similar -= 1
#                 similar_index.append(i)
#                 similar_index.append(j)
#     similar_index = list(set(similar_index))
#
#     for i in range(len(likelihood_parameter)):
#         if likelihood_parameter[i] <= ALPHA_THRESHOLD:
#             if i not in similar_index:
#                 count_similar -= 1
#
#     return count_similar


def write_SNPs_into_file(resDir, data, polyCounts_refLabel,
                         polyCounts_clustered,
                         proportion):
    res = defaultdict(list)
    nuclList = ['A', 'C', 'G', 'T']
    for item in data:
        tupleList = []

        pos = item[0]
        nucl = item[3:]
        nucl_ref = polyCounts_refLabel[pos]
        for i in range(len(nucl)):
            tupleList.append((nucl[i], nuclList[i]))
        tupleList.sort(key=itemgetter(0), reverse=True)

        # rewrite!!!
        SNP = []
        for one in tupleList:
            if one[1] == nucl_ref:
                continue
            else:
                if one[0] > polymorphic_cutoff:
                    SNP.append(one[1])

        for item in SNP:
            result_tmp = item + ',' + str(pos)
            nucl_label = polyCounts_clustered[pos]
            for label in nucl_label:
                res[label].append(result_tmp)

    proportion_dict = {}
    for item in proportion:
        proportion_dict.update({item[0]: item[2]})

    resName = resDir + '/haplotypes'
    with open('%s' % resName, 'w') as f:
        for key, value in res.items():
            # if key == 0:
            #     continue
            # else:
            key = str(key)
            proportion_item = proportion_dict[key]

            title = '>' + str(proportion_item)
            f.write('%s\n' % title)
            for item in value:
                f.write('%s\n' % item)


def generate_final_SNPs_single_hap(data, polyCounts_refLabel):
    res = defaultdict(list)
    nuclList = ['A', 'C', 'G', 'T']
    res = {'1.0':[]}
    for item in data:
        tupleList = []

        pos = item[0]
        nucl = item[3:]
        nucl_ref = polyCounts_refLabel[pos]
        for i in range(len(nucl)):
            tupleList.append((nucl[i], nuclList[i]))
        tupleList.sort(key=itemgetter(0), reverse=True)

        # rewrite!!!
        SNP = []
        for one in tupleList:
            if one[1] == nucl_ref:
                continue
            else:
                if one[0] > polymorphic_cutoff:
                    SNP.append(one[1])

        for item in SNP:
            result_tmp = item + ',' + str(pos)
            res['1.0'].append(result_tmp)

    return res


# filter position with total reads less than 10% of all position.
def filter(data):
    res = []

    # filter the position with total number of reads less than 10%
    total_read = [sum(item[1:]) for item in data]
    quantile = np.quantile(total_read, 0.1)
    total_read.sort()

    for item in data:
        if sum(item[1:]) <= quantile:
            continue
        else:
            res.append(item[:])

    return res


def run_EM(inFileName, numHaps, genomeFileLoc, resDir):
    # read polymorphic sites and filter
    data = read_polymorphicsSite(inFileName)
    count_before = len(data)
    data = filter(data)

    count_after = len(data)
    tmp = [count_before, count_after,
           round(count_after * 1.0 / count_before, 5)]
    print(inFileName + ' ratio of polymorphic sites survived after filter 10% of total reads: ' + ' '.join([str(item) for item in tmp]))

    # label polymorphic sites what is reference in that pos.
    polyCounts_refLabel = check_refernce_genome(data, genomeFileLoc)

    # normalized and filter individual letter
    data_normalized = normalized_filter(data, polyCounts_refLabel)

    # EM algorithm
    ROUND, pi, likelihood_parameter = EM(data_normalized, numHaps, resDir,
                                         polyCounts_refLabel)

    # # # hap=4
    # pi = [0.11245368519323012, 0.88754631480677]
    # likelihood_parameter = [[0.11245368519323012, 0.8875463148067699],
    #                         [0.88754631480677, 0.11245368519322996]]

    numHaps_predict = decide_num_haplotype(pi, likelihood_parameter)

    # #write EM parameter into files
    # resName = resDir + '/EM_parameter_pi'
    # with open('%s' %resName,'w') as f:
    #     tmp = '\n'.join([str(one) for one in pi])
    #     f.write('%s' %tmp)
    # resName = resDir + '/EM_parameter_alpha'
    # with open('%s' %resName,'w') as f:
    #     for item in likelihood_parameter:
    #         tmp = '\t'.join([str(one) for one in item])
    #         f.write('%s\n' %tmp)

    # get the posterior
    polyCounts_clustered, likelihood_sum = posterior_cluster(data_normalized,
                                                             pi,
                                                             likelihood_parameter,
                                                             numHaps,
                                                             polyCounts_refLabel)

    # compare the result with true label
    # compare_with_trueLabel(polyCounts_clustered,data,resDir)

    # calculate the proportion
    data_format = format_data(data, polyCounts_refLabel)
    proportion_normalized, proportion = calculate_proportion(data_format,
                                                             data_normalized,
                                                             polyCounts_refLabel,
                                                             polyCounts_clustered,
                                                             resDir)

    # write haplotypes SNPs into file
    write_SNPs_into_file(resDir, data_normalized, polyCounts_refLabel,
                         polyCounts_clustered, proportion_normalized)

    return numHaps_predict, likelihood_sum, len(data_normalized)


# automatically predict number of haplotypes
def automatica_predict_EM(parameter):
    print(parameter)
    inFileName = parameter[0]
    genomeFileLoc = parameter[1]
    resLoc = parameter[2]

    # print parameter
    previous_numHaps = -1
    likelihood_log = {}

    if not os.path.exists(resLoc):
        os.mkdir(resLoc)

    # calculate BIC when i=1
    BIC_value_1 = calculate_BIC_one_hap(inFileName, genomeFileLoc)

    # check whether it is only one haplotype
    flag_single_hap, snp_out_single = filter_one_hap(inFileName, genomeFileLoc)
    if flag_single_hap:
        print('No haplotype print')
        return 1, snp_out_single

    # predict 2 hap to 9 hap
    f1 = open(resLoc + '/strain_num_parameter.bed','w')
    for i in range(2, 10):
        # if i != 3:
        #     continue
        resDir = resLoc + '/hap_' + str(i)
        if not os.path.exists(resDir):
            os.mkdir(resDir)
        try:
            now_numHaps, likelihood_sum, nums_data = run_EM(inFileName, i,
                                                            genomeFileLoc,
                                                            resDir)
        except:
            print('##########################################################')
            print(parameter)
            print('#############################################################')
            break

        likelihood_log[i] = likelihood_sum
        flag_independent, bic_value = BIC(likelihood_log, i, nums_data)

        f1.write('strain_' + str(i) + '\t' + str(bic_value) + '\t' + str(BIC_value_1) + '\t' + str(now_numHaps) + '\t' + str(previous_numHaps) + '\t' + str(flag_independent) + '\n')
        if now_numHaps == 2:
            if bic_value >= BIC_value_1:
                print('No haplotype print')
                return 1, snp_out_single
        if now_numHaps == previous_numHaps or not flag_independent:
            break
        elif i == 2 and now_numHaps < i:
            return 2, snp_out_single
        else:
            previous_numHaps = now_numHaps

        # if now_numHaps == previous_numHaps or now_numHaps < i or not flag_independent:
        #     break
        # else:
        #     previous_numHaps = now_numHaps

    actual_numHaps = previous_numHaps
    print(actual_numHaps)
    f1.close()
    return actual_numHaps, snp_out_single

################################################################################
