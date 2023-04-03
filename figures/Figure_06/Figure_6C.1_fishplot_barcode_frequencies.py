#############################################################################################################################
##
##  RESCALING AND REFORMATING FOR FISHPLOT FOR D2C2R2
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################


import os
import math

path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6C/OUTPUT"
os.chdir(path)


fd = open("D2C2R2_CNV_subclones_programmable.txt","r")
fout = open("D2C2R2_CNV_subclones_programmable_rescale.csv","w")


line = fd.readline()
line = fd.readline()

### SET PARAMETERS
addon = 0.01  # ADJUSTMENT PARAMETER TO ACCOUNT FOR ROUNDING ERRORS
totsum = 100  # SCALE FACTOR
bottlesum = 20 # SCALE FACTOR FOR BOTTLENECK
roundnr = 5 # NR DIGITS WHEN ROUNDING

clones = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"] # ALL CLONES USED (INCL. CNV AND SUBCLONES)
times = ["0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"] # ALL TIME POINTS

timecloneDic = {}  # DICTIONARY CONTAINING EACH CLONE AND TIME POINT AND ASSOCIATED FREQUENCY

while line:
    line = line.split()
    time = line[3]
    clone = line[8]
    timeclone = time +"_"+ clone
    perc = float(line[6])
    perc1000 = perc*1000 # TAKE CLONE PERCENTAGE TIMES 1000
    try:
        log2 = math.log2(perc1000)  # SHIFT TO LOG2 SCALE TO BETTER VISUALIZE LOW FREQUENCY SUBCLONES
    except:
        log2 = 0

    log2 = round(log2, roundnr)
    if log2 < 0: # NEGATIVE VALUES ARE SET TO 0 (I.E. IF A CLONE HAS A VERY LOW FRACTION, IT WILL NOT BE VISUALIZED)
        log2 = 0
    timecloneDic[timeclone] = log2
    line = fd.readline()

### WRITE COLUMN WITH CLONE NUMBER
fout.write(",")
for c in clones:
    fout.write("C"+str(c)+",")
fout.write("\n")

###
# THE FISHPLOT R PACKAGE REQUIRES THAT PARENTAL CLONES HAVE AN EQUAL OR LARGER FRACTION THAN ITS CHILD CLONES
# THE CHILD CLONES HERE ARE MOSTLY THE ECB BASED SUBCLONES, AND THE PARENT CLONES THE CNV BASED CLONES
# DUE TO THE LOG2 RESCALING, WE NEED TO START COUNTING THE FREQUENCY OF EACH ECB SUBCLONE TO GET THE FREQUENCY OF THE CNV PARENT CLONES
# AND THEN SCALE ALL CLONES SO CLONE 1 - THE FOUNDING CLONES, HAS VALUE 100

for t in times:
    clone_five_sum = 0
    clone_four_sum = 0
    clone_two_sum = 0
    clone_one_sum = 0

    for c in clones:
        currTC = t+"_"+c
        try:
            currLog2 = timecloneDic[currTC]
        except:  #
            currLog2 = 0


        fives = ["11","12","13","14"]  # CNV CLONE 5 HAS CHILD CLONES 11,12,13,14
        fours = ["8","9","10"]  # CNV CLONE 4 HAS CHILD CLONES 8,9,10

        if c in fives:
            clone_five_sum += currLog2
        if c in fours:
            clone_four_sum += currLog2
        if c == "7":
            clone_two_sum += currLog2
        if c == "6":
            clone_one_sum += currLog2


    clone5 = clone_five_sum + addon  # CLONE FIVE FREQUENCY IS THE SUM OF THE CHILD CLONES PLUS "ADDON" TO ACCOUNT FOR ROUNDING ERRORS
    clone4 = clone5 + clone_four_sum + addon # CLONE FOUR FREQUENCY IS THE SUM OF THE CHILD CLONES INCL. CLONE 5
    clone3 = clone4 + addon # CLONE 3 HAS ONLY CLONE 4 AS A CHILD CLONE
    clone2 = clone3 + clone_two_sum + addon
    clone1 = clone2 + clone_one_sum + addon+0.00001 # AN EXTRA ADDON WAS REQUIRED HERE TO ACCOUNT FOR ROUNDING ERRORS

    scale = 100/clone1 # SCALE EVERYTHINGS SO CLONE 1 HAS VALUE 100

    if t == "3": # TO CREATE BOTTLE NECK VISUALIZATION OF ECB SELECTION, SCALE ON TIME POINT 3 WAS SET TO 20 FOR CLONE 1
        scale = 20 / clone1

    # WRITE CNV CLONES TO FILE
    fout.write(str(t)+",")

    fout.write(str(round(clone1*scale, roundnr))+","+str(round(clone2*scale, roundnr))+","+
               str(round(clone3*scale, roundnr))+","+str(round(clone4*scale, roundnr))+","+str(round(clone5*scale, roundnr)))

    # WRITE ECB SUBCLONES TO FILE
    for c in clones[5:]:
        currTC = t + "_" + c
        try:
            currLog2 = timecloneDic[currTC]
        except:
            currLog2 = 0
        currLog2_scale = round(currLog2*scale, roundnr)

        fout.write(","+str(currLog2_scale))

    fout.write(",")
    fout.write("\n")
