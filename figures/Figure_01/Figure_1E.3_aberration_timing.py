#############################################################################################################################
##
##  CALCULATE GENOME ABERRATION TIMING
##
##  OUTPUTS FOR EACH SAMPLE IF A CHROMOSOME ARM HAS A COPY NUMBER ABERRATION OR NOT
##  A CHROMOSOME ARM IS CALLED AS ABERRANT IF AT LEAST 25% OF WINDOWS HAVE A 25% INCREASE OR DECREASE COMAPARED TO MEDIAN NUMBER OF READS PER WINDOW
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

import os

### SET PATH
path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_1/"
os.chdir(path)

### GET ALL FILE NAMES
files = os.listdir(path+"INPUT/CNA_smooth_50kb_shallow_MA25_Centered_1/Summary")

### OPEN FILES TO WRITE IN
fout = open("PLOTS/Figure_1E_Timing_CNA/CNV_timing_per_sample.csv", "w")
fout.write("Sample"+","+"Donor"+","+"Clone"+","+"ECB_sample"+","+"Day"+","+"ChrArm"+","+"Direction"+","+"Aberration"+","+"Fraction_altered"+"\n")

### FUNCTION THAT CALLS IF A CHROMOSOME ARM IS ABERRANT OR NOT
def getAlts(file,name,fout):
    fd = open(path+"INPUT/CNA_smooth_50kb_shallow_MA25_Centered_1/Summary/"+file, "r")
    line = fd.readline()
    line = fd.readline()
    donor = name[0:2]
    clone = name[2:4]
    day = name[5:-1]
    ecb = "NO"

    while line:
        line = line.split()
        chrArm = line[0][3:]
        tot_win = float(line[3])
        win_alt_up = float(line[4])
        win_alt_down = float(line[5])

        if tot_win != 0:
            frac_up = win_alt_up/tot_win
            frac_down = win_alt_down/tot_win

            if float(frac_up) >= 0.25: # IF MORE THAN 25% OF WINDOWS HAVE INCREASED NUMBER OF READS THE CHROMOSOME ARM IS CALLED ABERRANT
                chrArmAlt = chrArm+"+"
                fout.write(name+","+str(donor)+","+str(clone)+","+ecb+","+str(day)+","+chrArm+","+"Up"+","+chrArmAlt+","+str(frac_up)+"\n")
            if float(frac_down) >= 0.25: # IF MORE THAN 25% OF WINDOWS HAVE DECREASED NUMBER OF READS THE CHROMOSOME ARM IS CALLED ABERRANT
                chrArmAlt = chrArm + "-"
                fout.write(name+","+str(donor)+","+str(clone)+","+ecb+","+str(day)+","+chrArm+","+"Down"+","+chrArmAlt+","+str(frac_down)+"\n")
        line = fd.readline()

### LOOP THROUGH ALL SAMPLES
for file in files:
    if "_50_moving" in file:
        if ".xlsx" not in file and "Mid" not in file:
            name,rest = file.split("_50_moving")
            if "WT" not in name and "Parent" not in name and "R" not in name: #Timing is not investigated for ECB cultures in figure 1
                getAlts(file,name,fout)
