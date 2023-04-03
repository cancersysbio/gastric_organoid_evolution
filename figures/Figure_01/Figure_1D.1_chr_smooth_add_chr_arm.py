#############################################################################################################################
##
##  ADD CHROMOSOME ARM LOCATION
##
##  AMENDS QDNASEQ OUTPUT FILES "CNV_SMOOT" WITH CHROMOSOME ARM LOCATION
##  USED FOR BOTH FIGURE 1 AND 2 (SHALLOW AND DEEP SEQUENCING)
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

import os

path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/"

### SET PARAMETERS
figure = "Figure_1" #Change here between Figure_1 and Figure_2
depth = "shallow" #Change here between deep and shallow
winsize ="50kb" #Change here between 1kb and 50kb

#figure = "Figure_2" #Change here between Figure_1 and Figure_2
#depth = "deep" #Change here between deep and shallow
#winsize ="1kb" #Change here between 1kb and 50kb

### GET ALL FILE NAMES
os.chdir(path)
files = os.listdir(path+figure+"/INPUT/CNA_smooth_"+winsize+"_"+depth) ### READ IN FILE NAMES

fd = open("sequtils/Supporting_files/hg38_cytoBand_simple_chr.txt", "r") ### READ IN GENOMIC POSITION FOR CHROMOSOME ARMS
line = fd.readline()
chrarmSplit = {}

### ADD CHROMOSOME CENTROMERE POSITION TO A DICTIONARY
while line:
    line = line.split()
    chr = line[3]
    arm = line[4]
    centro = int(line[2])
    if arm == "p":
        chrarmSplit[chr]=centro
    line = fd.readline()
fd.close()
print (chrarmSplit)

### ADD CENTROMERE POSITION AS LAST COLUMN TO THE "SMOOTH" OUTPUT FILE FROM QDNASEQ
def addChrAmr(file,name,chrarmSplit):
    fout = open(figure+"/INPUT/CNA_smooth_"+winsize+"_"+depth+"/"+name+"_CNA_smo_chrarm.txt", "w")
    fd = open(figure+"/INPUT/CNA_smooth_"+winsize+"_"+depth+"/"+file, "r")
    line = fd.readline()
    line = fd.readline()

    while line:
        line = line.split()
        for l in line:
            fout.write(str(l))
            fout.write("\t")
        chr = str(line[1])
        start = int(line[2])
        centro = chrarmSplit[chr]
        if start <= centro:
            chrarm = "p"
        if start > centro:
            chrarm = "q"
        fout.write(str(chr)+chrarm)
        fout.write("\n")
        line = fd.readline()


### LOOP THROUGH ALL "SMOOTH" FILES IN DIRECTORY
for file in files:
    if "CNA_smooth" in file:
        print(file)
        name,rest = file.split("_CNA_smooth")
        addChrAmr(file,name,chrarmSplit)


