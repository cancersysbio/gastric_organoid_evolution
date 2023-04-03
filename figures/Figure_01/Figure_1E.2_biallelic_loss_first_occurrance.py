#############################################################################################################################
##
##  BIALLELIC LOSS, TIMING FIRST OCCURRANCE
##
##  THE PROGRAM OUTPUTS THE FIRST OCCURANCE OF A BIALLELIC LOSS
##  ONLY ABERRATIONS THAT ARE PRESENT AT ALL SUBSEQUENT TIME POINTS ARE KEPT
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

import os

### SET PATH
path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_1/PLOTS/Figure_1E_Timing_CNA/"
os.chdir(path)

### OPEN FILES TO READ AND WRITE
fout = open("DualKO_summary_50kb_0.25_first_ocurrance_per_sample.csv","w")
fout.write("Donor" + "," + "Clone" + "," + "Day" + "," + "Aberration" + "," + "\n")

### FUNCTION TO IDENTIFY FIRST OCCURRANCE OF ABERRATION
def first_occurance(culture):
    fd = open("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/sequtils/Supporting_files/sWGS_samples_long_trajectory.csv", "r")
    line = fd.readline()
    line = fd.readline()
    days = []
    while line:
        line = line.split(",")
        donor = line[0]
        clone = line[1]
        day = line[2]
        DC = donor+clone
        if DC == culture:
            if int(day) not in days:
                days.append(int(day))
        line = fd.readline()
    fd.close()

    fd = open("DualKO_summary_50kb_0.25.csv", "r")
    line = fd.readline()
    line = fd.readline()
    abDic = {} ### STORES ALL DAYS AN ABERRATION OCCURRED

    ### LOOP THROUGH FILE WITH THAT CONTAINS TIMING OF ALL BIALLELIC LOSSES OF ALL SAMPLES
    while line:
        line = line.split(",")
        donor = line[1]
        clone = line[2]
        day = line[3]
        DC = donor+clone
        aberration = line[4].strip()
        if DC == culture: ### USE ONE CULTURE AT A TIME
            if aberration in abDic:
                currDays = abDic[aberration]
                updateDays = currDays + ";" + day
                abDic[aberration] = updateDays
            if aberration not in abDic:
                abDic[aberration] = day
        line = fd.readline()
    fd.close()
    abFinalDay = {}  ### STORE THE DAY OF FIRST OCCURRANCE OF THE ABERRATION
    droplist = []  ### IF ABERRATION NOT PRESENT AT ONE TIME POINT, STORE IN THIS LIST
    days.sort(reverse=True)
    ### LOOP THROUGH ALL DAYS IN REVERSED TEMPORAL ORDER
    for d in days:
        d = str(d)
        for aberration in abDic:
            abDays = abDic[aberration].split(";")
            if d in abDays:
                if aberration not in droplist:
                    abFinalDay[aberration] = d
            if d not in abDays:
                droplist.append(aberration)

    ### WRITE TO FILE
    for ab in abFinalDay:
        day = abFinalDay[ab]
        fout.write(culture[0:2]+","+culture[2:4]+","+str(day)+","+ab+","+"\n")

### RUN FOR ALL CULTURES
first_occurance("D1C1")
first_occurance("D1C2")
first_occurance("D1C3")
first_occurance("D2C1")
first_occurance("D2C2")
first_occurance("D2C3")
first_occurance("D3C1")
first_occurance("D3C2")
first_occurance("D3C3")

