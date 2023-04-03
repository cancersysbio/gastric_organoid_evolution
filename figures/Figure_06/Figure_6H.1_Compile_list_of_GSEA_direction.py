#############################################################################################################################
##
##  COMPILE LISTS OF GSEA OUTPUT VALUES FOR ECB REPLICATES
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################



import os
from os import listdir

path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6H/D2C2R2"
#path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6H/D2C1R1"
#path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6H/D2C1R2"
#path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6H/D3C2R1"


os.chdir(path)

fout = open("GSEA_compiled_all_direction_D2C2R2.txt","w")
#fout = open("GSEA_compiled_all_direction_D2C1R1.txt","w")
#fout = open("GSEA_compiled_all_direction_D2C1R2.txt","w")
#fout = open("GSEA_compiled_all_direction_D3C2R1.txt","w")

fout.write("Culture"+"\t"+"Timepoint"+"\t"+
           "Pathway"+"\t"+"Observed_Score"+"\t"+"p_val"+"\t"+"p_val_adj"+"\t"+"Leading_edge"+"\n")


for file in listdir(path):
    if "markers" in file:
        rest, sample, subclone, rest2, rest3  = file.split("_")
        if sample == "R2T2":
            sample = "D2C2R2"
            day = "173"
        if sample == "R1T2":
            sample = "D3C2R1"
            day = "189"
        if sample == "R1T7":
            sample = "D2C1R1"
            day = "258"
        if sample == "R2T7":
            sample = "D2C1R2"
            day = "258"
        if sample == "R1T12":
            sample = "D2C2R1"
            day = "315"
        if sample == "R2T12":
            sample = "D2C2R2"
            day = "315"
        if sample == "R3T12":
            sample = "D2C2R3"
            day = "315"

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

           fout.write(sample+"_"+subclone + "\t" + day + "\t")
           fout.write(pathway + "\t" +  str(value) + "\t" + direction + str(pval) + "\t" + direction + str(padj) + "\t"+ str(leadingEdge))
           fout.write("\n")
           line = fd.readline()






