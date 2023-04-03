#############################################################################################################################
##
##  INFERRED GROWTH TRAJECTORIES ECB CLONES LONG
##  Script to combine RGC matrix (barcode frequencies) and growth derivatives in long format for plotting
##
##
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

import os

### Set sample name and directory
sample_name = "D1C1" # Change sample name here
os.chdir("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/Growth_trajectory_ECB/Output")

### Open growth derivative and output files
fd = open(sample_name+"_Early_RG_abs_fitness_allReps.csv", "r")
fout = open("growth_deriv_heatmap_"+sample_name+"_Early_long_with_freq.txt", "w")

### Read in time points
line = fd.readline()
line = line.split(",")
tp = {}
growths = {}
n=0
for l in line:
    tp[n] = l.strip("\n")
    n+=1
line = fd.readline()

### Read in growth derivative per time point and RG into dictionary
while line:
    line = line.split(",")
    rg_rep = str(line[0][1:-1])
    n=1
    for l in line[1:]:
        time = str(tp[n][6:-1])
        rg_rep_tp = rg_rep + "_" + time
        l = l.strip("\n")
        if l == "":
            l=0
        growths[rg_rep_tp] = l
        n+=1
    line = fd.readline()
fd.close()

### Open RGC matrix file
matrixpath = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/Growth_trajectory_ECB/Input/norm_RGC_matrices/"
fd = open(matrixpath+sample_name+"_RGC_matrix_small_norm.txt", "r")

### Read in time points
line = fd.readline()
line = line.split("\t")
n=0
reptime = {}
for l in line:
    sample = l.split("_")
    go = "OK"
    rep = sample[1]
    if rep == "Parent":
        time = "T0"
        go = "NO"
    if rep == "Parent1":
        time = "T0"
        go = "NO"
    if rep == "Parent2":
        time = "T0"
        go = "NO"
    if go == "OK":

        time = sample[2].strip("\n")
    rt = rep+"_"+time
    reptime[n] = rt
    n+=1
line = fd.readline()

### Read in barcode frequency per time point and RG into dictionary
rg_reptime = {}
while line:
    line = line.split("\t")
    RG = line[0]
    n=0
    for l in line[1:]:
        rt = reptime[n]
        rgrt = str(RG)+"_"+str(rt)
        rg_reptime[rgrt] = l.strip("\n")
        n+=1

    line = fd.readline()

### Print all to file
fout.write("barcode"+"\t"+"RG"+"\t"+"rep"+"\t"+"RG_rep"+"\t"+"time"+"\t"+"Freq"+"\t"+"Deriv"+"\n")
for rg in growths:
    deriv = growths[rg]
    if rg in rg_reptime:
        freq = rg_reptime[rg]
    if rg not in rg_reptime:
        freq = "NA"
    rgs,rep,time = rg.split("_")
    rgs_rep = rgs+"_"+rep
    if freq != "Na":
        if deriv != "Na":
            fout.write(rg+"\t"+rgs+"\t"+rep+"\t"+rgs_rep+"\t"+time+"\t"+str(freq)+"\t"+str(deriv)+"\n")
