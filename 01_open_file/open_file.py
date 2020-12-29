

input_file_name = "codon_transfer.txt"
#output_file_name = "test.txt"

# #https://ithelp.ithome.com.tw/articles/10161708

file = open(input_file_name, mode = "r" )
# output_file = open(output_file_name , mode = "w")  

#### Reading each lines into list
#print( file.readlines() )
import re
for i in file.readlines():
    readline = i.rstrip()
    #print(readline)
    if re.findall(r'(AAA)\t(.+)',readline):
        
        print(re.findall(r'AAA',readline))
    ## output_file.write(readline+"\n")

file.close()
# output_file.close()



i = ""
readline = ""

# with open(input_file_name, mode = "r") as file:
#     for i in file:
#         readline = i.rstrip()   ## rstrip() 去換行
#         print(readline)



