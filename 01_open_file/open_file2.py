import sys
import re

input_file_name = sys.argv[1]
i = ""
readline = ""
count = 0 

readline_list = []
amino_acid = {}
codon_table = {}

with open(input_file_name,mode = "r") as file:
    for i in file:
        readline = i.rstrip()
        count = count+1
        readline_list = readline.split('\t')
        #print(readline_list)
        if len(readline_list) < 3:
            pass
        else:
            
            #print(readline_list[0])
            try:
                codon_table[readline_list[0]] = (readline_list[2])
            except KeyError:
                codon_table[readline_list[0]] = (readline_list[2])
            
            try:
                amino_acid[readline_list[1]] = int(amino_acid[readline_list[1]])+1
            except KeyError:
                amino_acid[readline_list[1]] = 1
                #pass
            
            #print(amino_acid[readline_list[1]])
        #print(readline)

# for i in amino_acid.keys():
#     print(str(i) + "\t" + str(amino_acid[i]))

# print("Total: " + str(count) + " codon")
# print("Total: " + str(len(amino_acid.keys()) -1 ) + " amino acid & stop codon" )


input_file_name = sys.argv[2]
with open(input_file_name,mode = "r") as file:
    for i in file:
        readline = i.rstrip()
        
        if re.findall( r"^>", readline):
            print(readline)
        else:
            pass
            seq = readline
            #print(readline)
            protein_seq = ""
            for i in range(0, len(seq), 3): 
                codon = seq[i:i + 3] 
                try:
                    protein_seq+=codon_table[codon]
                except KeyError:
                    protein_seq+=codon_table[codon]
                sep=""
                #print((codon_table[codon][0]))
                #protein = sep.join(protein_seq)
            print((protein_seq))
        # #count = count+1
        #readline_list = readline.split('\t')











