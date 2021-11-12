import argparse

from get_sample_reads_frequency import get_reads_frequency
from merge_filter_polymorphic_sites import getMergedPolymorphicSites
from normalized_filter_polymorphic_sites import normalized_matrix
from mixture_model import *
from get_final_strain import *
from n1_nodes_filter import *

parser = argparse.ArgumentParser()
# group1 = parser.add_argument_group('General')
parser.add_argument('--output_name', help='Give a unique output name',
                    default=None)
parser.add_argument('--genome_len', help='Input genome length', default=None)
parser.add_argument('--genome_name', help='Input genome name', default=None)
parser.add_argument('--genome_file_loc', help='Input genome file location', default=None)
parser.add_argument('--bam_loc_file', help='Input bam file loc list', default=None)
parser.add_argument('--res_dir', help='result directory', default=None)
#
args = parser.parse_args()

#result directory
if not args.res_dir:
    script_path = os.path.dirname(os.path.realpath(__file__))
else:
    script_path = os.path.dirname(args.res_dir)

genome_len = int(args.genome_len)
genome_name = args.genome_name
bam_loc_file = args.bam_loc_file
res_dir = os.path.abspath(args.res_dir) + '/'
genome_file = os.path.abspath(args.genome_file_loc)
sample_name = args.output_name

res_tmp_dir = res_dir + sample_name
if not os.path.exists(res_tmp_dir):
    command = 'mkdir -p ' + res_tmp_dir
    os.system(command)

# get sample reads frequency ####
if not os.path.exists(res_tmp_dir + '/sample_reads_frequency'):
    command = 'mkdir -p ' + res_tmp_dir + '/sample_reads_frequency'
    os.system(command)
get_reads_frequency(genome_len,genome_name,bam_loc_file,res_tmp_dir + '/sample_reads_frequency')

# get merged polymorphic sites
print('************** Getting merged filter polymorphic sites*************')
reads_frequency_file_list = [res_tmp_dir + '/sample_reads_frequency/' + i for i in os.listdir(res_tmp_dir + '/sample_reads_frequency')]
getMergedPolymorphicSites(genome_len,reads_frequency_file_list, res_tmp_dir,genome_file)

# normalized merged polymorphic sites
print('************** normalize merged polymorhpic sites')
normalized_matrix(res_tmp_dir + '/merged_filter_polymorphic_sites',genome_file + '_updated')

######## remove the n1 nodes that sum of vector smaller than 5 ###############

merged_file = res_tmp_dir + '/merged_filter_polymorphic_sites_normalized_100'
frequency_file_loc = res_tmp_dir + '/sample_reads_frequency'
result_dir = res_tmp_dir
get_n1_nodes_filter(merged_file,frequency_file_loc,result_dir)

# running mixture to get R strain ########################################
run_mixture(res_tmp_dir,genome_file,sample_name)




### get final strain ##############
merged_file = res_tmp_dir + '/merged_filter_polymorphic_sites_normalized_100_filter'
frequency_file_loc = res_tmp_dir + '/sample_reads_frequency_filter'
result_dir = res_tmp_dir + '/result'
mixture_strain_file = res_tmp_dir + '/' + sample_name + '_haplotypes_0.6'
if os.path.exists(result_dir) == False:
    os.makedirs(result_dir)

import time
#
start = time.time()
strain_dict = get_strain(merged_file, mixture_strain_file, frequency_file_loc,result_dir)
end = time.time()
print('time used')
print(end - start)

time_file = result_dir + '/time_used.bed'

f1 = open(time_file,'w')

f1.write(str(end-start) + '\n')

f1.close()


