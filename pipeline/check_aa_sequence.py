##############################################################
# 
# checking variation postition whather change amino acid sequence
# 
# usage : python3 check_aa_sequence.pl $ARGV[1]  $ARGV[2] 
# 
# python3 check_aa_sequence.pl ATCC824_reorganized.fna mutation_in_coding_region.txt 
# 
# 
# 	$ARGV[1] >>> reference cds position file (***_cds_from_genomic.fna)
# 	$ARGV[2] >>> input vcf file (mutation_in_coding_region.txt )
# 
# 

import sys
import re


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



genome={}

###############################
####                       ####
####  Reading genome fasta ####
####                       ####
###############################

# >NC_003030.1 Clostridium acetobutylicum ATCC 824 chromosome, complete genome
# TCAAGAAAAACGATTCTGCCAACCTTTAGCAATGAGATTTTATTCCTTACATATTAAATTTTTCACTTCTGTTGATAAAT
file_path= sys.argv[1]
with open(file_path, mode = "r", encoding = "utf8") as file:
    for i in file:
        
        readline = i.rstrip()  ## rstrip() >>> 去除換行符號
        
        
        if re.match( "^>", readline):
            # pass
            read_element = readline.split('\t')
            Chromosome_id = read_element[0]
            Chromosome_id =Chromosome_id.replace(">","")
            try:
                genome[Chromosome_id] = ""
            except KeyError:
                genome[Chromosome_id] = ""
            
        else:
            genome[Chromosome_id] = (readline)


protein ="" 


###############################
####                       ####
####  Reading codon table  ####
####                       ####
###############################
def translate(seq): 
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        
        'ata':'I', 'atc':'I', 'att':'I', 'atg':'M', 
        'aca':'T', 'acc':'T', 'acg':'T', 'act':'T', 
        'aac':'N', 'aat':'N', 'aaa':'K', 'aag':'K', 
        'agc':'S', 'agt':'S', 'aga':'r', 'agg':'R', 
        'cta':'L', 'ctc':'L', 'ctg':'L', 'ctt':'L', 
        'cca':'P', 'ccc':'P', 'ccg':'P', 'cct':'P', 
        'cac':'H', 'cat':'H', 'caa':'Q', 'cag':'Q', 
        'cga':'R', 'cgc':'R', 'cgg':'R', 'cgt':'R', 
        'gta':'V', 'gtc':'V', 'gtg':'V', 'gtt':'V', 
        'gca':'A', 'gcc':'A', 'gcg':'A', 'gct':'A', 
        'gac':'D', 'gat':'D', 'gaa':'E', 'gag':'E', 
        'gga':'G', 'ggc':'G', 'ggg':'G', 'ggt':'G', 
        'tca':'S', 'tcc':'S', 'tcg':'S', 'tct':'S', 
        'ttc':'F', 'ttt':'F', 'tta':'L', 'ttg':'L', 
        'tac':'Y', 'tat':'Y', 'taa':'_', 'tag':'_', 
        'tgc':'C', 'tgt':'C', 'tga':'_', 'tgg':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein += table[codon] 
            #print(protein,"")
    return protein 




# Chromosome_id	mutation_pos	alt	ref	mutation_type	cds_start	cds_end	transcription_way	locus_taq	old_locus_taq	protein_id	description
# NC_003030.1	196943	TAAAAA	TAAAAAA	deletion	195829	196990	+	CA_RS00990	CA_C0173	WP_013913506.1	sensor histidine kinase
# NC_003030.1	303842	G	T	snp	303812	305674	+	CA_RS01510	CA_C0272	WP_010963595.1	APC family permease
# NC_003030.1	345995	T	C	snp	345987	346103	+	CA_RS01670	CA_Cr027		5S ribosomal RNA
# NC_003030.1	346085	G	C	snp	345987	346103	+	CA_RS01670	CA_Cr027		5S ribosomal RNA
# NC_003030.1	389758	A	T	snp	389645	390223	+	CA_RS01900	CA_C0331	WP_010963654.1	ECF transporter S component
# NC_003030.1	391277	G	A	snp	390304	391371	-	CA_RS01905	CA_C0332	WP_013913668.1	glycoside hydrolase family 26 protein
# NC_003030.1	423814	GAAAAAA	GAAAAAAA	deletion	423328	424036	+	CA_RS02070	CA_C0364	WP_013913511.1	MBL fold metallo-hydrolase
# NC_003030.1	423881	G	A	snp	423328	424036	+	CA_RS02070	CA_C0364	WP_013913511.1	MBL fold metallo-hydrolase




#############################
####                     ####
####  Reading input file ####
####                     ####
#############################
vcf_path= sys.argv[3]
with open(vcf_path, mode = "r", encoding = "utf8") as file:
    for i in file:
        
        readline = i.rstrip()  ## rstrip() >>> 去除換行符號
        
        read_element = readline.split('\t')
        
        Chromosome_id = read_element[0]
        mutation_pos = read_element[1]
        alt = read_element[2]
        ref = read_element[3]
        mutation_type = read_element[4]
        cds_start = read_element[5]
        cds_end = read_element[6]
        transcription_way = read_element[7]
        locus_taq = read_element[8]
        old_locus_taq = read_element[9]
        protein_id = read_element[10]
        description = read_element[11]
        output = read_element[12]
        
        
        
        
        #





