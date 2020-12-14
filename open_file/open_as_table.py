import sys
import pandas as pd


input_file_name = sys.argv[1]

i = ""
readline = ""

sheet = pd.read_table(input_file_name, header = None, sep = '\t' )
file_len = sheet.index.stop
column_list = sheet.columns

# print(sheet[0])
# print(sheet.index)
# print(column_list)

for i in list(range(sheet.index.stop)):
    DNA_codon = sheet[0][i]
    protein = sheet[1][i]
    print(DNA_codon+"\t"+protein)
    #print()


