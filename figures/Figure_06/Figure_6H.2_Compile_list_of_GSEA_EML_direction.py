#############################################################################################################################
##
##  COMPILE LISTS OF GSEA OUTPUT VALUES FOR EML SAMPLES
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################


import os
from os import listdir

path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6H/EML"

os.chdir(path)

fout = open("GSEA_compiled_all_direction_EML.txt","w")

fout.write("Culture"+"\t"+"Timepoint"+"\t"+
           "Pathway"+"\t"+"Observed_Score"+"\t"+"p_val"+"\t"+"p_val_adj"+"\t"+"Leading_edge"+"\n")


for file in listdir(path):

    if "markers" in file:
        print (file)
        rest, sample, timepoint, rest2, rest3, rest4, rest5  = file.split("_")

        fd = open(file,"r")
        line = fd.readline()
        line = fd.readline()

        while line:
           line = line.split()
           pathway = line[0]
           value = float(line[1])
           pval = float(line[2])
           padj = float(line[3])
           leadingEdge = line[4]

           if value < 0:
               direction = "-"
           if value > 0:
               direction = ""

           if pval == 0:
               pval =  0.00000000000001
           if padj == 0:
               padj = 0.00000000000001

           print (str(value))

           fout.write(sample + "\t" + timepoint+ "\t")
           fout.write(pathway + "\t" + str(value) + "\t" + direction + str(pval) + "\t" + direction + str(padj) + "\t"+ str(leadingEdge))
           fout.write("\n")
           line = fd.readline()






