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
            read_element = readline.split(' ')
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
        
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
        'tac':'Y', 'tat':'Y', 'taa':'*', 'tag':'*', 
        'tgc':'C', 'tgt':'C', 'tga':'*', 'tgg':'W', 
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


mutation_info = {}      ## key = (chr_id+","+pos), value = (ref+","+alt)
cds_info = {}           ## key = "cds id", value = Chromosome_id+","+cds_start+","+cds_end+","+transcription_way+","+locus_taq+","+old_locus_taq+","+protein_id+","+description+","+output"
cds_all_mutation = {}   ## key = "cds id", value = (pos1, pos2, ...)




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

mut_key = ""
#

#############################
####                     ####
####  Reading input file ####
####                     ####
#############################
vcf_path= sys.argv[2]
with open(vcf_path, mode = "r", encoding = "utf8") as file:
    for i in file:
        
        readline = i.rstrip()  ## rstrip() >>> 去除換行符號
        
        read_element = readline.split('\t')
        if str(read_element[0]) == "Chromosome_id":
            pass
        else:
            Chromosome_id = str(read_element[0])
            mutation_pos = str(read_element[1])
            alt = str(read_element[2])
            ref = str(read_element[3])
            mutation_type = str(read_element[4])
            cds_start = str(read_element[5])
            cds_end = str(read_element[6])
            transcription_way = str(read_element[7])
            locus_taq = str(read_element[8])
            
            try:
                old_locus_taq = str(read_element[9])
            except IndexError:
                old_locus_taq = ""
            
            try:
                protein_id = str(read_element[10])
            except IndexError:
                protein_id = ""
            
            description = str(read_element[11])
            
            mut_key=(Chromosome_id+","+mutation_pos)
            
            
            #########################
            ####                 ####
            #### 建立dictionary  ####
            ####                 ####
            #########################
            try:
                mutation_info[mut_key] = (alt+"\t"+ref)     
            except KeyError:
                mutation_info[mut_key] = (alt+"\t"+ref)
            
            try:
                cds_info[locus_taq] = str(Chromosome_id+"\t"+cds_start+"\t"+cds_end+"\t"+transcription_way+"\t"+locus_taq+"\t"+old_locus_taq+"\t"+protein_id+"\t"+description)
            except KeyError:
                cds_info[locus_taq] = str(Chromosome_id+"\t"+cds_start+"\t"+cds_end+"\t"+transcription_way+"\t"+locus_taq+"\t"+old_locus_taq+"\t"+protein_id+"\t"+description)
            
            try:
                cds_all_mutation[locus_taq] = (cds_all_mutation[locus_taq]+"\t"+str(mut_key))
            except KeyError:
                cds_all_mutation[locus_taq] = (str(mut_key))
            #


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])



########################################
####                                ####
#### 最後整理確認cds 中所有mutation  ####
####                                ####
########################################

info = ""
fragment1 = ""
fragment2 = ""
fragment3 = ""


for i in list( cds_all_mutation.keys() ):
    
    
    

    info = str(cds_info[i])
    
    read_element = info.split('\t')
    Chromosome_id = read_element[0]
    cds_start = (int(read_element[1])-1)
    cds_end = (int(read_element[2]))
    
    
    # print(cds_start)
    
    transcription_way = read_element[3]
    locus_taq = read_element[4]
    old_locus_taq = read_element[5]
    protein_id = read_element[6]
    description = read_element[7]
    
    
    genome_fna = str(genome[Chromosome_id])
    
    
    ###################################
    ####                          #####
    ####  Reference DNA sequence  #####
    ####                          #####
    ###################################
    cds_fna = genome_fna[cds_start:cds_end]
    mut_cds_fna = cds_fna
    
    # print(">"+locus_taq+"\t"+old_locus_taq+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+description+"\t"+transcription_way+"\n"+cds_fna)
    # cds_all_mutation[i] = cds_all_mutation[i].replace(",", "", 1)
    mutation_list = []
    mutation_list = cds_all_mutation[i].split("\t")
    # print (mutation_list)
    indel_base = 0
    
    for mut_key in mutation_list:
        # print (mut_key)
        
        
        ref_alt = mutation_info[mut_key].split("\t")
        # print (ref_alt)
        ref = ref_alt[1]
        alt = ref_alt[0]
        
        pos_list = mut_key.split(",")
        pos_1 = ( (int(float(pos_list[1]) - cds_start -1) + indel_base))
        pos_2 = (pos_1+len(ref))
        
        if len(ref) == len(alt):    ## SNP
            pos_2 = (pos_1+len(ref) )
        elif len(ref) > len(alt):   ## deletion
            pass
        elif len(ref) < len(alt):   ## insertion
            pass
        
        
        
        indel_base = (indel_base + (len(ref) - len(alt)))
        
        fragment1 = mut_cds_fna[0:pos_1]
        fragment2 = mut_cds_fna[pos_1:pos_2]
        fragment3 = mut_cds_fna[pos_2:]
        
        
        
        
        if fragment2 == ref:
            fragment2 = alt
        else:
            pass
            print(mut_key+"\tSomethong wrong!!!"+"\t"+str(fragment2)+"\t"+str(ref))
        
        
        mut_cds_fna = (fragment1+fragment2+fragment3)
    
    if transcription_way == "-":
        cds_fna = reverse_complement(cds_fna)
        mut_cds_fna = reverse_complement(mut_cds_fna)
    
    # print(">"+locus_taq+"\t"+old_locus_taq+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+description+"\t"+transcription_way)
    # print(mut_key)
    # print(cds_fna+"\n>alt\n"+mut_cds_fna)

    #

