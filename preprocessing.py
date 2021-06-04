import argparse
import os, sys


'''
python preprocessing.py --sample_name test --pair1 ./test_data/test.read1.fastq --pair2 ./tes
t_data/test.read2.fastq --process 6 --genome_name ./test_data/GCF_000016525.1_ASM1652v1_genomic.fna --res_dir ./test_res_data
'''


# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument('--sample_name', help='Give a unique sample name',
                    default=None)
parser.add_argument('--pair1', help='Give a unique sample name',
                    default=True)
parser.add_argument('--pair2', help='Give a unique sample name',
                    default=True)
parser.add_argument('--process', help='Give a unique sample name',
                    default=True)
parser.add_argument('--genome_name', help='Input genome name', default=None)
parser.add_argument('--res_dir', help='result directory', default=None)

args = parser.parse_args()

script_dir = os.path.dirname(os.path.abspath(__file__))

#result directory
res_dir = os.path.abspath(args.res_dir)



if not os.path.exists(args.res_dir):
    command = 'mkdir -p ' + res_dir
    os.system(command)


#Trim the fastq
res_forward_read1 = res_dir + '/' + args.sample_name + \
                    '_forward.filter.fastq'
res_reverse_read2 = res_dir + '/' + args.sample_name + \
                    '_reverse.filter.fastq'
res_forward_unpair_read1 = res_dir + args.sample_name + \
                           '_forward_unpair.fastq'
res_reverse_unpair_read2 = res_dir + args.sample_name + \
                           '_reverse_pair.fastq'
trimmomatic_path = script_dir + '/trimmomatic-0.36.jar'
command1 = ' '.join(['java', '-jar', trimmomatic_path, 'PE', '-phred33',
                     args.pair1, args.pair2, res_forward_read1,
                     res_forward_unpair_read1, res_reverse_read2,
                     res_reverse_unpair_read2,
                     'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10', 'LEADING:3',
                     'TRAILING:3', 'SLIDINGWINDOW:4:20', 'MINLEN:36'])
print('###########Step1#######################################################')
print(command1)
os.system(command1)

#################### build index for ref genome ###################

os.system('bowtie2-build %s %s' % (args.genome_name,args.genome_name))

####################################################################

#map to sam
res_sam = res_dir + '/' + args.sample_name + '.sam'
command2 = ' '.join(['bowtie2', '-p', args.process, '-N', '1', '-X',
                     '1000', '-x', args.genome_name, '-1', res_forward_read1,
                     '-2', res_reverse_read2, '-S', res_sam])
print('###########Step2#######################################################')
print(command2)
os.system(command2)


#sam to bam
res_bam = res_dir + '/' + args.sample_name + '.bam'
command3 = ' '.join(['samtools', 'view', '-Sb', res_sam, '>', res_bam])
print('###########Step3#######################################################')
print(command3)
os.system(command3)

#sort bam
res_sort_bam = res_dir + '/' + args.sample_name + '_sorted.bam'
command4 = ' '.join(['samtools', 'sort', '-o', res_sort_bam, res_bam])
print('###########Step4#######################################################')
print(command4)
os.system(command4)


#bam index
command5 = ' '.join(['samtools', 'index', res_sort_bam])
print('###########Step5#######################################################')
print(command5)
os.system(command5)


os.remove(res_bam)
os.remove(res_sam)
os.remove(res_forward_read1)
os.remove(res_reverse_read2)
os.remove(res_forward_unpair_read1)
os.remove(res_reverse_unpair_read2)

