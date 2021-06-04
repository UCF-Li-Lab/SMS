
import os.path
import Constants
from EM import automatica_predict_EM
import sys

############################################################################
SEQUENCE_ERROR = Constants.SEQUENCE_ERROR
############################################################################

#EM algorithm
def run_mixture(res_tmp_dir,genome_file,sample_name):
    res_tmp_dir = res_tmp_dir + '/'
    polymorphic_file_name = res_tmp_dir + 'merged_filter_polymorphic_sites_100_filter'
    em_dir = res_tmp_dir + 'res_em/'
    num_haps_predict, snp_output_single = automatica_predict_EM((
        polymorphic_file_name, genome_file, em_dir))

    if num_haps_predict == 1 or num_haps_predict == -1:
        command = 'rm -r ' + em_dir
        os.system(command)

        #write haplotype for single into file
        resName = res_tmp_dir + sample_name + '_haplotypes'
        with open('%s' % resName, 'w') as f:
            for key, value in snp_output_single.items():
                title = '>' + str(1.0)
                f.write('%s\n' % title)
                for item in value:
                    f.write('%s\n' % item)

        sys.exit('no haplotypes found')


    res_hap = em_dir + 'hap_' + str(num_haps_predict) + '/haplotypes'
    command = 'cp ' + res_hap + ' ' + res_tmp_dir + sample_name + '_haplotypes_0.6'
    os.system(command)


    #remove intermediate files
    #command = 'rm -r ' + em_dir
    #os.system(command)

    print('mixture running Done')


# res_tmp_dir = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/xin_195_data/testing_result/NC_000962.3_test'
# genome_file = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/xin_195_data/original_genome/GCF_000195955.2_ASM19595v2_genomic.fna'
# sample_name = 'NC_000962.3'
#
# run_mixture(res_tmp_dir,genome_file,sample_name)