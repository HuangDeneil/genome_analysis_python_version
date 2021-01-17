###################################
# 
# pyhton3 genome_file_reorganized.py $ARGV[0]
# 
# 
# 
# 
# 
# ##

import sys
import re



file_path= sys.argv[1]

genome_out = open("reorganized_genome.txt", mode = "w" )

count = 0
with open(file_path, mode = "r", encoding = "utf8") as file:
    for i in file:
        
        readline = i.rstrip()  ## rstrip() >>> 去除換行符號
        if re.match( "^>", readline):
            
            if count == 0 :
                genome_out.write(readline+"\n")
            else:    
                genome_out.write("\n"+readline+"\n")
            
        else:
            genome_out.write(readline)
        
        count = 1




genome_out.close()
