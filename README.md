# SMS

SMS: a novel approach for microbial strain genome reconstruction in multiple samples


=====================================================================================

SMS: a novel approach for microbial strain genome reconstruction in multiple samples

Author: Saidi Wang

Date: 06/04/2021

=====================================================================================


                                 ################# Working OS #################
                                 
Our model is working on Linux Operating System.

                                 ################# Virtual Environment #################
                                 
Most Easy way to set up environment is using conda. 

A MSMS_py36.yml file was provided in the source directory.

(1) Install Conda:

https://docs.conda.io/projects/conda/en/latest/user-guide/install/

(2) Open a terminal in mixture model directory, create virtual environment by command:

conda env create -f ./MSMS_py36.yml

(3) Before running model, go into the virtual environment by command:

conda activate MSMS_py36

                               ###### Common problem with Virtual Environment #########
                               
After you install some package, if our model are running with error. It is mostly caused by the rpy2 package, you probably manually install R, and then install rpy2.

(1) Inside the environment created above, Conda install R

conda install -c r r

(2) Inside the environment created above, install rpy2 with specific version

pip install rpy2

                              ###### Manual Install Without Virtual Environment #######
                              
If you do not like conda or can not use it by conda for some reason, then in your own environment, the following prerequisite software with required version should be installed:

(1) python>=3.6.10 If there is no python  installed, you can download and install python from (http://www.python.org/download/). You can use "python -V" command to check whether python is installed and the version of Python.

(2) Install R

https://www.r-project.org/ or install it by conda(https://anaconda.org/r/r)

(3) pysam=0.15.3

Detail to install pysam in the link: https://pysam.readthedocs.io/en/latest/installation.html. Or use three

commands below one by one:

conda config --add channels r

conda config --add channels bioconda

conda install pysam

(4) rpy2==2.9.4

conda install rpy2

(5) jenkspy=0.1.5

pip install jenkspy

(6) scipy=1.4.1

conda install scipy=1.4.1

(7) Bowtie2 2.3 or newer (optional if you are only using our preprocessing script)

Conda install bowtie2

(8) Samtools 1.9 (optional if you are using our preprocessing script)

Conda install samtools=1.9 or https://anaconda.org/bioconda/samtools

=====================================================================================

                                  ###### Preprocessing #######
                                  
Since MSMS require bam format input, you had better input the sorted Bam file as input. However, if you only have FastQ format data, you can get a hint from our common preprocessing pipeline. But please check your experimental specification before applying our preprocessing step, because our preprocessing may be not suitable for your specific experiments. So you may need to adjust our preprocessing script based on your own experiments or trim read protocol, like single-end read, trim customized bar code or others.

(1) Use conda virtual environment in Virtual Environment:

(2) Manually prerequisite

If you would like to use our preprocessing steps, there are two more additional prerequisite, see (7) and (8) in Manual Install Without Virtual Environment.

(3) Sample command:

python preprocessing.py --sample_name test --pair1 ./example_test_data/test.read1.fastq -- pair2 ./example_test_data/test.read2.fastq --process 6 --
genome_name ./example_test_data/ref.fna --res_dir ./test_res_data

--sample_name: sample name

--pair1: forward read for FastQ

--pair2: reverse read for FastQ

--process: # of processor to run

--genome_name: reference genome

--res_dir : result directory

=====================================================================================

                          ###### Parameters for Run MSMS model #######
(1) Go into virtual environment if you deploy by virtual environment

conda activate MSMS_py36

(2) Run our MSMS model:

python ./MSMS.py

--output_name output result name: Given a unique sample name for running result

--genome_len genomeLength: the corresponding genome length

--genome_name genomeName: the genome name(name inside FASTA file)

--genome_file_loc genomeLoc: genome FASTA file location. Absolute file location is better

--bam_loc_file bam_file_location: the files contain all the sorted bam file which have been mapped to the reference genome.

--res_dir directory_name: the result directory name

(3) sample input

The test data has been saved in the directory of ‘example’, so you can just navigate into directory where MSMS.py located. And use below command to run example 

data. The result will save into directory ‘example’.

python MSMS.py  --output_name group1_1_NC_009515.1 --genome_len 1853160 --genome_name NC_009515.1 --genome_file_loc example/reference/GCF_000016525.1_ASM1652v1_genomic.fna --bam_loc_file  example/bam_result/group1_1_NC_009515.1_bam_loc.bed --res_dir example/MSMS_result

=====================================================================================

                                      ###### Interpret results #######
                                      
1.output_name/sample_reads_frequency_filter directory: it’s the filter location and its frequencies in each sample.

2. output_name/result

‘MSMS_result.bed’ --------

This is the final prediction of strains of MSMS model. It includes coverage with corresponding polymorphic sites. In the example below, there are two strains predicted. Each strain starts with '>' and its corresponding coverage. From next row, they are the polymorphic sites of above strain. Each row is represented by allele and its position. In the example below, two strains with coverage 39.95 and 10.29. 

Example result (part):

>strain_39.95

T,1399289

A,918579

G,202077

A,1030566

C,1731682

G,1701635

C,641623

T,1585394

T,1195555

T,110682

T,1540799

>strain_10.29

A,672356

T,1394505

A,1588695

A,1037717

C,534803

A,1162319

A,169157

T,1485071

A,1120883

T,145207

G,56444

T,1532716

A,305246

C,340455

‘strain_in_each_sample.txt’  -------------- it’s the information of the strains occurred in each sample.

‘1’ ---- the strain occurs in the sample.

‘0’ ----- the strain did not occur in the sample

example result:

sample_name	 strain_0_39.95	strain_1_10.29	strain_2_28.35	strain_3_19.71

NC_009515.1_1a2_0.01_100_400_0.86a19.47_0.001_2264a8463_sorted.bam_reads_frequency	0	1	0	1

NC_009515.1_1a3_0.01_100_400_2.19a30.0_0.001_9528a3289_sorted.bam_reads_frequency	0	1	1	0

NC_009515.1_1a4_0.01_100_400_3.33a21.16_0.001_719a6069_sorted.bam_reads_frequency	1	1	1	0

NC_009515.1_1a4_0.01_100_400_3.62a6.16_0.001_8402a3220_sorted.bam_reads_frequency	1	1	0	0

NC_009515.1_2a4_0.01_100_400_0.53a12.68_0.001_5387a1108_sorted.bam_reads_frequency	1	0	0	1

‘time_used.bed’ ------------- it’s the time used to get the strains from the bam file( the final SMS method step).

========================================================================

                  ####### Contact Information ######

Please do not hesitate to reach out to me if you have questions.

Saidi Wang (tjwangsaidi@knights.ucf.edu)
Haiyan Nancy Hu (haihu@cs.ucf.edu)
Xiaoman Shawn Li (xiaoman@mail.ucf.edu)
