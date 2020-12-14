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
snp_dict = {}
insertion_dict = {}
deletion_dict = {}

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
                    snp_dict[key]=info
                except KeyError:
                    pass
            elif len(alt) > len(ref):
                try:
                    insertion_dict[key]=info
                except KeyError:
                    pass
            elif len(alt) < len(ref):
                try:
                    deletion_dict[key]=info
                except KeyError:
                    pass
           
            # try :
            #     print (deletion_dict[key])
            # except KeyError:
            #     pass




'''
with open (open_path ,mode = "r", encoding = "utf8" ) as file:
    for i in file:
        tradition_cn = sc2tc(i).rstrip()           ## rstrip() >>>  perl chomp;
        #print(tradition_cn )
        #pattern = '樣本型別：'

    regex = re.compile(r'[Vv][Ii][Rr][Uu][Ss]')  # if match virus
    m = regex.search(top_type)



'''


gene_locus=[]
tmp_list=[]
locus_tag_chrmosome_dict={}
locus_tag_info_dict={}
protein_id=""
protein_id_dict={}

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
            position=(read_element[3]+"-"+read_element[4])
            if re.findall( r"gene", read_element[2]):
                
                info_match = re.findall( r"\;old_locus_tag\=(.+)", read_element[8])
                if info_match:
                    gene_locus = info_match[0].split(';')
                    # if re.findall( r";", str(info_match[0])):
                    #     print (gene_locus[0])
                    try:
                        locus_tag_chrmosome_dict[position] = (chrmosome+"-"+gene_locus[0])
                    except KeyError:
                        locus_tag_chrmosome_dict[position] = (chrmosome+"-"+gene_locus[0])
                
                #
            elif re.findall( r"(CDS|RNA)", read_element[2]):
                
                ## grep protein id
                info_match = re.findall( r"\;protein_id=(.+)?", read_element[8])
                if info_match:
                    tmp_list = info_match[0].split(';')
                    protein_id = tmp_list[0]
                    try:
                        protein_id_dict[position] = iprotein_id
                    except KeyError:
                        protein_id_dict[position] = iprotein_id
                
                
                # 
                if re.findall( r"RNA", read_element[2]):
                    try:
                        non_coding_gene[position] = True
                    except KeyError:
                        non_coding_gene[position] = True
                
                
                # pass





