

## Project Introduction
This project mainly aims to standardize the editing information uploaded by users, in order to accurately locate the position of the editing area on the genome.According to the type of gene sequence file uploaded by users, this system mainly supports the parsing and operation of two types of files.

1.When a user uploads a genome file as an fna file and a matching editing information csv file, the system will perform a blast alignment based on the first 100 bp sequence of the user specified editing area to accurately locate it on the genome. Finally, the standard output results after localization will be output to serve downstream work
    For downstream scenarios: it will be divided into i.designing only sgRNA (only_sgRNA) ii.designing both sgRNA and related primers (both_sgRNA_primer)
    i:  data1 = {
                "input_file_path":"./input/editor_info.csv",
                "ref_genome":"./input/GCA_000011325.1_ASM1132v1_genomic.fna",
                "data_preprocessing_workdir":"/home/XXX/tmp/data_preprocessing/output/",
                "scene":"only_sgRNA"
            }
    ii: data2 = {
                "input_file_path":"./input/editor_info123.csv",
                "ref_genome":"./input/GCA_000011325.1_ASM1132v1_genomic.fna",
                "data_preprocessing_workdir":"/home/XXX/tmp/data_preprocessing/output/",
                "scene":"both_sgRNA_primer"
            }
    iii: data3 = {
                "input_file_path":"./input/sgRNA_editing_info.csv",
                "ref_genome":"./input/GCA_000011325.1_ASM1132v1_genomic.fna",
                "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                "scene":"only_primer",  
            }

    output: '/home/XXX/tmp/data_preprocessing/output//info_input.csv'

2.When a user uploads a genome gb file and a matching editing information csv file, the system will first convert it into a genome fna file, then extract the editing region sequence based on the editing region provided by the user, and finally output the standardization results to serve downstream work.
    For downstream scenarios: it will be divided into i.designing only sgRNA (only_sgRNA) ii.designing both sgRNA and related primers (both_sgRNA_primer)

    i:  data4 = {
                "input_file_path":"./input/4-21-input.csv",
                "ref_genome":"./input/eco.gb",
                "data_preprocessing_workdir":"/home/XXX/tmp/data_preprocessing/output/",
                "scene":"only_sgRNA",
            }
    ii: data5 = {
                "input_file_path":"./input/4-20-input.csv",
                "ref_genome":"./input/eco.gb",
                "data_preprocessing_workdir":"/home/XXX/tmp/data_preprocessing/output/",
                "scene":"both_sgRNA_primer",
            }
    iii: data6 = {
                "input_file_path":"./input/4-23-input.csv",
                "ref_genome":"./input/eco.gb",
                "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                "scene":"only_primer",  
            }
            
    output: ['/home/yanghe/tmp/data_preprocessing/output/info_input.csv', '/home/yanghe/tmp/data_preprocessing/output/eco.fna']

## product enviroment
 python 3.8
## Software
$ pip install requirements.txt

## use
python parse_input_to_df