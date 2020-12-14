

input_file_name = "codon_transfer.txt"
i = ""
readline = ""

output_file = open("test.txt", mode = "w")  

with open(input_file_name, mode = "r") as file:
    for i in file:
        readline = i.rstrip()   ## rstrip() 去換行
        #print(readline)
        output_file.write(readline+"\n")

output_file.close()

# #https://ithelp.ithome.com.tw/articles/10161708
# file = open(input_file_name, mode = "r" )

# ## readlines in list
# #print( file.readlines() )

# for i in file.readlines():
#     readline = i.rstrip()
#     print(readline)

# file.close()
