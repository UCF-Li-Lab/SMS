
import pysam
import os

def fetch_sites_reads(genomeLen, genomeName, in_bam_file, res_dir,sample_name):
    # calculate the ACGT frequency for each position
    # initial
    # genomeLen = 1853160
    # genomeName = 'NC_009515.1'
    # inName = '/media/student/study_working/project/project12/point1/data/simulated_data/simulatedReads_sortedBam/default/0_sorted.bam'

    print('generate polymorphic sites')
    genomeLen = int(genomeLen)
    res = []
    for _ in range(genomeLen):
        res.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})


    samfile = pysam.AlignmentFile('%s' % in_bam_file, 'rb')
    posID = 0
    # genome start location 0, the following pileup parameter is the same.
    for pileupColumn in samfile.pileup(genomeName, 0, genomeLen):
    #for pileupColumn in samfile.pileup(genomeName, 1000000, 1001000):
        if posID % 100000 == 0:
            print(posID)
        posID += 1

        pos = pileupColumn.pos
        count = pileupColumn.n
        tmp = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for pileupread in pileupColumn.pileups:
            if (not pileupread.is_del) and (not pileupread.is_refskip):
                nucl = pileupread.alignment.query_sequence[pileupread.query_position]
                nucl = nucl.upper()
                tmp[nucl] += 1

        for key, value in tmp.items():
            res[pos][key] = value

    samfile.close()

    # calculate polymorphics sites
    nuclOrder = ['A', 'C', 'G', 'T']
    polymorphicSites = []
    for loc in range(len(res)):
        pos_nucl = res[loc]
        formatOutput = [loc]

        # flag_nonZero = 0
        for nucl in nuclOrder:
            nucl_count = pos_nucl[nucl]
            formatOutput.append(nucl_count)
            # if nucl_count != 0:
            #     flag_nonZero += 1

        # if flag_nonZero >= 1:
        polymorphicSites.append(formatOutput)

    res_file_name = res_dir + '/' + sample_name + '_reads_frequency'
    with open('%s' %res_file_name, 'w') as f:
        for item in polymorphicSites:
            one = '\t'.join([str(one) for one in item])
            f.write('%s\n' %one)
    return polymorphicSites


def get_reads_frequency(genomeLen,genomeName,in_bam_loc_file,res_dir):
    ### getting reads count #####
    print('###########Step1 get sample reads frequency ##########################')
    f1 = open(in_bam_loc_file, 'r')
    lines1 = f1.readlines()
    bam_file_list = []
    for line1 in lines1:
        bam_file_list.append(line1.strip())  ### bam_file

    for in_bam_file in bam_file_list:
        sample_name = os.path.split(in_bam_file)[-1]
        fetch_sites_reads(genomeLen, genomeName, in_bam_file, res_dir, sample_name)