#############################################################################################################################
##
##  INFERRED GROWTH TRAJECTORIES LONG
##  Script to combine growth derivatives, passaging numbers and growth fold change in long format for plotting
##
##
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################


import os
from os import listdir

### Set path and time point parameter
path1 = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_3/Figure_3B_Growth_trajectory_nonECB/"
path2 = path1+"Final_Output/passage_num/" # folder with passage corrections
os.chdir(path1)
timepoint = "Late" # Mid or Late

### Open growth derivative and output files
fd = open("Final_Output/deriv_by_samp_EarlyTo"+timepoint+".csv","r")
fout = open("Final_Output/nonECB_growth_deriv_heatmap_Early"+timepoint+"_raw_counts.csv", "w")

### Read in growth derivative per sample and passage into a dictionary
line = fd.readline()
line = line.split(",")
samples = {}
growths = {}
n=0
for l in line:
    samples[n] = l.strip("\n")
    n+=1

line = fd.readline()
while line:
    line = line.split(",")
    #print (line)
    passage = line[0]
    n=1
    for l in line[1:]:
        sample = samples[n]
        samplepass = str(sample) + ":" + passage
        l = l.strip("\n")
        if l == "":
            l=0
        growths[samplepass] = l
        n+=1
    line = fd.readline()
fd.close()

### Read in cell counts and fold change growth (in cell number between passages)
fd = open("Final_Input/Organoid_logbook_"+timepoint+"_Trajectory_Adjusted.csv", "r") # LATE
#fd = open("Final_Input/Organoid_logbook_"+timepoint+"_Trajectory.csv", "r") # MID

line = fd.readline()
line = fd.readline()

counts = {}

while line:
    line = line.split(",")
    sample_pass = "D"+line[0]+"C"+line[1] + ":" + line[2]
    count = line[4]
    growth_FC = line[7].strip("\n")
    count_FC = count+"_"+growth_FC
    counts[sample_pass] = count_FC
    line = fd.readline()

### Read in passage corrections (first couple of passages after TP53 KO were not counted - and are here ignored)
passage_correction = {}

for file in listdir(path2):
    if "_Early"+timepoint+"_" in file:
        temp = file.split("_Early"+timepoint+"_")
        sample = temp[1][:-4]
        fd = open(path2 + file,"r")
        line = fd.readline()
        line = fd.readline()
        while line:
            line = line.split(",")
            sample_pass = sample + ":" + line[0]
            passage_correction[sample_pass] = line[1].strip("\n")
            line = fd.readline()

### Print all to file in long format
fout.write("sample,oriPass,newPass,growthDeriv,count,growth_FC"+"\n")

for c in counts:
    count,growth_FC = counts[c].split("_")
    sample,realPass = c.split(":")

    if c in passage_correction:  # c=sample + passage
        newPass = passage_correction[c]
        newSample = sample +":" + newPass
        growthDeriv = growths[newSample]
        fout.write(str(sample)+",")
        fout.write(str(realPass)+",")
        fout.write(str(newPass)+",")
        fout.write(str(growthDeriv)+",")
        fout.write(str(count)+",")
        fout.write(str(growth_FC)+"\n")
