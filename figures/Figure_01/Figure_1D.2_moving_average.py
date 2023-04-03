#############################################################################################################################
##
##  MOVING AVERAGE
##
##  ADD MOVING AVERAGE TO QDNASEQ OUTPUT FILES AND TRANSFORMS READS PER WINDOW TO READ COUNTS CENTERED AROUND 1 (INSTEAD OF LOG 2 VALUE)
##  USED FOR BOTH FIGURE 1 AND 2 (SHALLOW AND DEEP SEQUENCING)
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

import os
import numpy as np

### SET PARAMETERS
figure = "Figure_1" #Change here between Figure_1 and Figure_2
depth = "shallow" #Change here between deep and shallow
winsize ="50kb" #Change here between 1kb and 50kb

#figure = "Figure_2" #Change here between Figure_1 and Figure_2
#depth = "deep" #Change here between deep and shallow
#winsize ="1kb" #Change here between 1kb and 50kb


### GET ALL FILE NAMES
path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/"
os.chdir(path)
files = os.listdir(path+figure+"/INPUT/CNA_smooth_"+winsize+"_"+depth) ### READ IN FILE NAMES

### CHECK IF PATH FOR NEW FILES EXISTS, OTHERWISE CREATE

# Check whether the specified path exists or not
isExist = os.path.exists(figure + "/INPUT/CNA_smooth_" + winsize + "_" + depth + "_MA25_Centered_1")
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(figure + "/INPUT/CNA_smooth_" + winsize + "_" + depth + "_MA25_Centered_1")

isExist = os.path.exists(figure + "/INPUT/CNA_smooth_" + winsize + "_" + depth + "_MA25_Centered_1/Summary")
if not isExist:
    os.makedirs(figure + "/INPUT/CNA_smooth_" + winsize + "_" + depth + "_MA25_Centered_1/Summary")

### FUNCTION FOR MOVING AVERAGE
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

### TEST MOVING AVERAGE FUNCTION
test = [0,1,2,3,4,5,6,7,8,7,6,5,5,5]
print (running_mean(test,5))

### FUNCTION THAT OUTPUTS
def getAve(file,name):
    # FILE THAT CALCULATES NUMBER OF ALTERED WINDOWS PER FIL
    # A WINDOW IS CALLED AS ALTERED IF IT HAS A 25% INCREASE OR DECREASE COMAPARED TO MEDIAN NUMBER OF READS PER WINDOW

    fout = open(figure+"/INPUT/CNA_smooth_"+winsize+"_"+depth+"_MA25_Centered_1/Summary/"+name+"_moving_25bp_average_summary_windows_altered_25percent.txt", "w")
    # FILE WITH MOVING AVERAGE ADJUSTED WINDOWS
    fout2 = open(figure+"/INPUT/CNA_smooth_"+winsize+"_"+depth+"_MA25_Centered_1/"+name + "_moving_25bp_average.txt", "w")
    fd = open("sequtils/Supporting_files/hg38_cytoBand_simple_chr.txt", "r") ### READ IN GENOMIC POSITION FOR CHROMOSOME ARMS

    line = fd.readline()
    fout.write("ChrArm" + "\t" + "Start" + "\t" + "End" + "\t" + "total_win" + "\t" + "win_alt_up" + "\t" + "win_alt_down" + "\n")
    fout2.write("Region" + "\t" + "Chr" + "\t" + "Start" + "\t" + "End" + "\t" + "ChrArm"+ "\t" + "log2_val" + "\t" +
        "Norm_reads_center_1" + "\n")

    while line: # LOOP THROUGH CYTOBAND FILE
        line = line.split()
        chr = line[0][:-1]
        start = int(line[1])
        end = int(line[2])
        fd2 = open(figure+"/INPUT/CNA_smooth_"+winsize+"_"+depth+"/"+file, "r")
        line2 = fd2.readline()
        line2 = fd2.readline()
        inloop = "FALSE"
        nr_win = 0
        nr_win_up = 0
        nr_win_down = 0
        chrarm_vals = []
        gencoord = []
        #print chr
        while line2: # LOOP THROUGH "SMOOTH" QDNASEQ OUTPUT FILE
            line2=line2.split()
            l2_coord = line2[0]
            l2_chr = line2[1]
            l2_start = line2[2]
            l2_end = line2[3]
            val = float(line2[4])
            try:
                l2_chrarm = line2[5] # SKIP WINDOWS THAT DON'T HAVE CHROMOSOME ARM ANNOTATION
            except:
                print (line2)
            l2_comb = l2_coord+"_"+l2_chr+"_"+l2_start+"_"+l2_end+"_"+str(val)+"_"+l2_chrarm

            if chr[3:] == l2_chr:
                if start <= int(l2_start) and end > int(l2_start):
                    val2 = 2**val # SHIFT FROM LOG 2 VALUE TO COPY NUMBER VALUE
                    gencoord.append(l2_comb) # STORE GENOMIC COORDINATES
                    chrarm_vals.append(val2) # STORE COPY NUMBER VALUES
                    inloop = "TRUE"
            if chr[3:] != l2_chr and inloop == "TRUE": # BREAK LOOP IF CHROMOSOME DIFFERS BETWEEN THE CYTOBAND AND QDNASEQ OUTPUT FILES
                break
            line2 = fd2.readline()
        fd2.close()
        if inloop == "TRUE":
            chrarm_vals_runMean = running_mean(chrarm_vals,25) # CALCULATE RUNNING AVERAGE ON 25 WINDOWS
            gencoord_short = gencoord[12:-12] # REMOVE FIRST AND LAST 12 WINDOWS PER CHROMOSOME

            for n in range(len(gencoord_short)):
                # PRINT RUNNING AVERAGE VALUE PER WINDOW
                this_chord,this_chr,this_start,this_end,this_val,this_chrarm = gencoord[n].split("_")
                running_mean_val = chrarm_vals_runMean[n]
                fout2.write(this_chord+"\t"+this_chr+"\t"+this_start+"\t"+this_end+"\t"+this_chrarm+"\t"+this_val+"\t"+
                            str(running_mean_val)+"\n")

            # CALCULATE NUMBER OF WINDOWS ALTERED PER CULTURE
            for vals in chrarm_vals_runMean:
                nr_win += 1
                if vals >= 1.25:  # 25 percent increase
                    nr_win_up += 1
                if vals <= 0.75:  # 25 percent decrease
                    nr_win_down += 1

            # PRINT NUMBER OF WINDOWS ALTERED PER CULTURE TO FILE
            for l in line[0:3]:
                fout.write(str(l))
                fout.write("\t")
            fout.write(str(nr_win))
            fout.write("\t")
            fout.write(str(nr_win_up))
            fout.write("\t")
            fout.write(str(nr_win_down))
            fout.write("\n")
        line = fd.readline()
    fd.close()
    fout.close()

# LOOP THROUGH ALL "SMOOTH" QDNASEQ OUTPUT FILES
for file in files:
    if "_CNA_smo_chrarm" in file:
        print(file)
        name,rest = file.split("_CNA_smo_chrarm")
        getAve(file,name)

