##############################################################
# 
# checking variation postition whather at protein coding regein
# 
# usage : perl check_snp_position.py $ARGV[0] $ARGV[1] 
# 
# perl check_snp_position.py GCF_000008765.1_ASM876v1_genomic.gff HOL_S10_L001.bwa.sorted.calls.filted_005.vcf 
# 
# 
# 	$ARGV[0] >>> reference cds position file (***_cds_from_genomic.fna)
# 	$ARGV[1] >>> input vcf file
# 
# 
# 


import sys
import re

vcf_path = ""  ## input file
readline = ""  ## file reading in line
read_element = [] ## each line element

chro = ""
pos = ""
alt = ""
ref = ""

key = ""
info = ""

chro_dict= {}
type_dict= {}
non_coding_gene= {}
###############################
###                         ###
###    read vcf position    ###
###                         ###
###############################
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HOL1.bwa.sorted.bam

vcf_path= sys.argv[1]
with open(vcf_path, mode = "r", encoding = "utf8") as file:
    for i in file:
        
        readline = i.rstrip()  ## rstrip() >>> 去除換行符號
        
        if re.match( "^#", readline):
                readline = "header"
                #print (re.match( "^#", readline).group(0))
        else:
            read_element = readline.split('\t')
            chro = read_element[0]
            pos = read_element[1]
            ref = read_element[3]
            alt = read_element[4]
           
            key = (chro+"\t"+pos)
            info = (chro+"\t"+pos+"\t"+alt+"\t"+ref)
           
            if len(alt) == len(ref):
                try:
                    chro_dict[key]=info
                    type_dict[key]=("snp")
                except KeyError:
                    chro_dict[key]=info
                    type_dict[key]=("snp")
            elif len(alt) > len(ref):
                try:
                    chro_dict[key]=info
                    type_dict[key]=("insertion")
                except KeyError:
                    chro_dict[key]=info
                    type_dict[key]=("insertion")
            elif len(alt) < len(ref):
                try:
                    chro_dict[key]=info
                    type_dict[key]=("deletion")
                except KeyError:
                    chro_dict[key]=info
                    type_dict[key]=("deletion")
            



read_element = []
tmp_list=[]

chrmosome = ""
position_start = ""
position_end = ""
direction = ""
gene_locus = ""
old_gene_locus = ""
description = ""
protein_id = ""

position = "" ## `核心` >>> 所有dictionary 共用 key (cds region)

transcription_way = {}
locus_tag_dict = {}
old_locus_tag_dict = {}
protein_id_dict = {}
description_dict = {}

## 放最後突變位子發生地點
mutation_location = {} 


###########################
###                     ###
###    read gff file    ###
###                     ###
###########################
gff_path= sys.argv[2]
with open(gff_path, mode = "r", encoding = "utf8") as file:
    for i in file:
        
        readline = i.rstrip()  ## rstrip() >>> 去除換行符號
        if re.match( r"^#", readline):
            readline = "header"
            #print (re.match( "^#", readline).group(0))
        else:
            read_element = readline.split('\t')
            chrmosome = read_element[0]
            
            position_start = read_element[3]
            position_end = read_element[4]
            position = (chrmosome+"\t"+position_start+"\t"+position_end)
            
            if read_element[2] == "region":
                readline = "header"
            else:
                try:
                    transcription_way[position] = read_element[6]
                except KeyError:
                    transcription_way[position] = read_element[6]
            
            if re.findall( r"gene", read_element[2]):
                #######################
                #####             #####
                #####  抓gene id  #####
                #####             #####
                #######################
                info_match = re.findall(r"\;locus_tag=(.+)", read_element[8])
                if info_match:
                    tmp_list = info_match[0].split(';')
                    try:
                        locus_tag_dict[position] = (tmp_list[0])
                    except KeyError:
                        locus_tag_dict[position] = (tmp_list[0])
                
                ###########################
                #####                 #####
                #####  抓old gene id  #####
                #####                 #####
                ###########################
                info_match = re.findall( r"\;old_locus_tag\=(.+)", read_element[8])
                if info_match:
                    tmp_list = info_match[0].split(';')
                    try:
                        old_locus_tag_dict[position] = (tmp_list[0])
                    except KeyError:
                        old_locus_tag_dict[position] = (tmp_list[0])
                
                
            elif re.findall( r"(CDS|RNA)", read_element[2]):
                
                ##  protein id
                info_match = re.findall( r"RefSeq\:(.+)?", read_element[8])
                if info_match:
                    tmp_list = info_match[0].split(';')
                    protein_id = tmp_list[0]
                    try:
                        protein_id_dict[position] = protein_id
                    except KeyError:
                        protein_id_dict[position] = protein_id
                
                ##  protein description
                info_match = re.findall(r"\;product=(.+)?", read_element[8])
                if info_match:
                    tmp_list = info_match[0].split(';')
                    description = tmp_list[0]
                    try:
                        description_dict[position] = description
                    except KeyError:
                        description_dict[position] = description
                
                if re.findall( r"RNA", read_element[2]):
                    try:
                        non_coding_gene[position] = True
                    except KeyError:
                        non_coding_gene[position] = True
        
        #######################################
        #####                             #####
        #####    確認突變位置是否在 cds 中  #####
        #####                             #####
        #######################################
        if readline == "header" :
            pass
        else:
            for k in chro_dict.keys():
                read_element = k.split('\t')
                chro = read_element[0]
                pos = read_element[1]
                
                if chro == chrmosome:
                    #position = (chrmosome+"\t"+position_start+"\t"+position_end)
                    if int(position_start) <= int(pos) :
                        if int(pos) <= int(position_end) :
                            try:
                                mutation_location[k] = position
                            except KeyError:
                                mutation_location[k] = position
                            #print (position+"\t"+k)





########################
#####              #####
##### 清空暫存變數  #####
#####              #####
########################
key = ""
info = ""
mutation_type = ""
position = ""
position_start = ""
position_end = ""
direction = ""
gene_locus = ""
old_gene_locus = ""
protein_id = ""
description = ""
read_element = []



###########################
#####                 #####
#####  Output result  #####
#####                 #####
###########################
print ("Chromosome_id\tmutation_pos\talt\tref\tmutation_type\tcds_start\tcds_end\ttranscription_way\tlocus_taq\told_locus_taq\tprotein_id\tdescription")
for k in chro_dict.keys():
    try:
        position = mutation_location[k]
    except KeyError:
        position = ""
    
    key = k
    try:
        info = chro_dict[key]
    except KeyError:
        info = ""
        
    try:
        mutation_type = type_dict[key]
    except KeyError:
        mutation_type = ""

    
    ###################################
    ####                           ####
    ####   加入cds information     ####
    ####                           ####
    ###################################
    if position == "": ### 空值代表 突變不發生在 transcription 區域
        #pass
        position_start = ""
        position_end = ""
        direction = ""
        gene_locus = ""
        old_gene_locus = ""
        protein_id = ""
        description = ""
    else:
        #####################################
        ####                             ####
        ####  非空值代表是coding region   ####
        ####                             ####
        #####################################
        read_element = position.split('\t')
        position_start = read_element[1]
        position_end = read_element[2]
        
        
        
        try:
            direction = transcription_way[position]
        except KeyError:
            direction = ""
        
        try:
            gene_locus = locus_tag_dict[position]
        except KeyError:
            gene_locus = ""
        
        try:
            old_gene_locus = old_locus_tag_dict[position]
        except KeyError:
            old_gene_locus = ""
        
        try:
            protein_id = protein_id_dict[position]
        except KeyError:
            protein_id = ""
        
        try:
            description = description_dict[position]
        except KeyError:
            description = ""
        
    print(info+"\t"+mutation_type+"\t"+position_start+"\t"+position_end+"\t"+direction+"\t"+gene_locus+"\t"+old_gene_locus+"\t"+protein_id+"\t"+description)
        
        #






