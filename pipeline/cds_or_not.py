#############################################################
# This program is checking :
#   variation postition whather in protein coding regein or not.
# 
# usage : python3 cds_or_not.py $ARGV[1] 
# 
#  	$ARGV[1] >>> cds_position_primary_check.txt
# 
# python3 cds_or_not.py cds_position_primary_check.txt
# 
# 
# 

import sys
import re


# Chromosome_id	mutation_pos	alt	ref	mutation_type	cds_start	cds_end	transcription_way	locus_taq	old_locus_taq	protein_id	description
# NC_003030.1	196943	TAAAAA	TAAAAAA	deletion	195829	196990	+	CA_RS00990	CA_C0173	WP_013913506.1	sensor histidine kinase
# NC_003030.1	303842	G	T	snp	303812	305674	+	CA_RS01510	CA_C0272	WP_010963595.1	APC family permease
# NC_003030.1	310671	T	TA	deletion							
# NC_003030.1	345919	G	A	snp							
# NC_003030.1	345995	T	C	snp	345987	346103	+	CA_RS01670	CA_Cr027		5S ribosomal RNA
# NC_003030.1	346085	G	C	snp	345987	346103	+	CA_RS01670	CA_Cr027		5S ribosomal RNA
# NC_003030.1	389758	A	T	snp	389645	390223	+	CA_RS01900	CA_C0331	WP_010963654.1	ECF transporter S component
# NC_003030.1	391277	G	A	snp	390304	391371	-	CA_RS01905	CA_C0332	WP_013913668.1	glycoside hydrolase family 26 protein

Chromosome_id = ""
mutation_pos = ""
alt = ""
ref = ""
mutation_type = ""
cds_start = ""
cds_end = ""
transcription_way = ""
locus_taq = ""
old_locus_taq = ""
protein_id = ""
description = ""
output = ""


file_cds = open("in_cds.txt", mode = "w" )

readline = ""
read_element = []

vcf_path= sys.argv[1]
with open(vcf_path, mode = "r", encoding = "utf8") as file:
    for i in file:
        
        readline = i.rstrip()  ## rstrip() >>> 去除換行符號
        
        read_element = []
        read_element = readline.split('\t')
        
        Chromosome_id = read_element[0]
        mutation_pos = read_element[1]
        alt = read_element[2]
        ref = read_element[3]
        mutation_type = read_element[4]
        
        try:
            # cds_start = str(read_element[5])
            locus_taq = read_element[8]
        except IndexError:
            # cds_start = ""
            locus_taq = ""
        
        
        if locus_taq == "" :
            pass
        else:
            output = readline
            file_cds.write (output + "\n" )
        
        
        #




file_cds.close()












