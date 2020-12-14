import sys



input_file_name = sys.argv[1]
i = ""
readline = ""
count = 0 

readline_list = []
amino_acid = {}

with open(input_file_name,mode = "r") as file:
    for i in file:
        readline = i.rstrip()
        count = count+1
        readline_list = readline.split('\t')
        #print(readline_list)
        if len(readline_list) < 3:
            pass
        else:
            try:
                amino_acid[readline_list[1]] = int(amino_acid[readline_list[1]])+1
            except KeyError:
                amino_acid[readline_list[1]] = 1
                #pass
            
            #print(amino_acid[readline_list[1]])
        #print(readline)

for i in amino_acid.keys():
    print(str(i) + "\t" + str(amino_acid[i]))

print("Total: " + str(count) + " codon")
print("Total: " + str(len(amino_acid.keys()) -1 ) + " amino acid & stop codon" )




