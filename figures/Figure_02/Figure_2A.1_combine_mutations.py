#############################################################################################################################
##
##  COMBINE MUTATIONS
##
##  Combines mutations from all .maf files
##
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################


import os
from os import listdir

path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/Mutect2_multiallele"

os.chdir(path)

fout = open("All_mutations_combined.txt","w")

for file in listdir(path):
    n = 0
    if "Mutect2_multiallele_pass.vcf.maf" in file:
        if "R2T" not in file:
            sample,temp = file.split("_Mutect2_")
            print(file, sample)

            fd = open(file, "r")
            line = fd.readline()
            line = fd.readline()

            if n == 1:
                line = fd.readline()
            n = 1

            while line:
                line = line.split("\t")
                fout.write(sample)
                for l in line[0:13]:
                    fout.write("\t" + str(l))
                for l in line[39:45]:
                    fout.write("\t" + str(l))
                fout.write("\n")
                line = fd.readline()
            fd.close()