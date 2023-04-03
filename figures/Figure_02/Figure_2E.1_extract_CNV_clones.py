#############################################################################################################################
##
##  EXTRACT COPY NUMBER ALTERATION SUBCLONES
##
##  THE PROGRAM INPUTS MEDIAN CENTERED AND MOVING AVERAGED SMOOTHENED QDNASEQ OUTPUT FILES
##  AND OUTPUTS REGIONS THAT ARE AT LEAST 100 WINDOWS LONG AND WHERE EACH WINDOW HAS AT LEAST A 12.5%
##  READ COUNT CHANGE COMPARED TO MEDIAN READ COUNT
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

import os
#path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_1/INPUT/CNA_smooth_50kb_shallow_MA25_Centered_1"
#path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/INPUT/CNA_smooth_50kb_deep_MA25_Centered_1"
path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/INPUT/CNA_smooth_50kb_shallow_MA25_Centered_1"

os.chdir(path)
files = os.listdir(path)

### SET PARAMETERS
limitUp = 1.125 # Specify percentage upregulated to count as a CNA clone
limitDown = 0.875 # Specify percentage downregulated to count as a CNA clone
minLenRegion = 100 # Specify minimal region length

### LOOP THROUGH EACH SAMPLE
def getVariable(sample):
    fout = open("CNV_clones_12.5p_q100/"+sample + "_CNV_clones.txt", "w") # OUTPUT FILE
    fout.write("window"+"\t"+"sample"+"\t"+"chr"+"\t"+"start"+"\t"+"end"+"\t"+"nr_win"+"\t"+"aveCopyNo"+"\n")
    fd = open(sample + "_50_moving_25bp_average.txt", "r") # READ IN MOVING AVERAGE QDNASEQ FILE
    line = fd.readline()
    line = fd.readline()
    n=0
    altChr = "XXX" # KEEPS TRACK OF IF THE PREVIOUS LINE IS PART OF THE SAME CHROMOSOME OR NOT
    altEnd = 0 # KEEPS TRACK OF THE LAST WINDOW TO PRINT IF E.G. DIRECTION CHANGES
    altDirection = "none" # KEEPS TRACK OF ALTERATION DIRECTION (UP OR DOWN)
    q = 0 # NUMBER OF ALTERED WINDOWS IN A ROW
    vals = 0 # AVERAGE ALTERATION VALUE OF ALL WINDOWS IN A CNA BASED SUBCLONE
    while line:
        line = line.split()
        chr = str(line[1])
        start = line[2]
        end = line[3]
        val = float(line[6])

        # GET DIRECTION OF ALTERATION (IF ALTERATION GREATER THAN THRESHHOLD)
        if val >= limitUp:
            direction = "up"
        if val <= limitDown:
            direction = "down"
        if val > limitDown and val < limitUp:
            direction = "none"

        # IF THE ALTERATION FOR A WINDOW IS LARGER THAN TRESHHOLD AND THE DIRECTION IS THE SAME AS LAST WINDOW
        # THEN KEEP TRACK OF NUMBER OF WINDOWS ALTERED
        if direction == "up" or direction == "down":  #
            if altChr == chr and direction == altDirection:
                if q == 1:
                    storeStart = start
                q += 1
                vals += val
                altEnd = end

        # PRINT ONLY IF THERE HAS BEEN A CHANGE (E.G. NEW CHROMOSOME OR DIRECTION CHANGED FROM UP TO NONE)
        printa = "FALSE"
        if altChr != chr:
            printa = "TRUE"
        if val > limitDown and val < limitUp:
            printa = "TRUE"
        if direction != altDirection:
            printa = "TRUE"

        # PRINT OUT CNA CLONE ONLY IF THE NUMBER OF WINDOWS ALTERED IS LARGER THAN 100
        if printa == "TRUE":
            if q >= minLenRegion:
                if vals / float(q) > limitUp or vals / float(q) < limitDown:
                    fout.write(str(n))
                    fout.write("\t")
                    fout.write(str(sample))
                    fout.write("\t")
                    fout.write(str(altChr))
                    fout.write("\t")
                    fout.write(str(storeStart))
                    fout.write("\t")
                    fout.write(str(altEnd))
                    fout.write("\t")
                    fout.write(str(q))
                    fout.write("\t")
                    fout.write(str(vals/float(q)))
                    fout.write("\n")
                    q = 1
                    vals = val

        if val > limitDown and val < limitUp:
            q = 1
            vals = val

        altStart = start
        altEnd = end
        altChr = chr
        altDirection = direction
        n+=1
        line = fd.readline()
    fd.close()

    # PRINT THE LAST REGION TOO IF THE NUMBER OF WINDOWS ARE MORE THAN 100
    if q >= minLenRegion:
        fout.write(str(n))
        fout.write("\t")
        fout.write(str(sample))
        fout.write("\t")
        fout.write(str(altChr))
        fout.write("\t")
        fout.write(str(storeStart))
        fout.write("\t")
        fout.write(str(altEnd))
        fout.write("\t")
        fout.write(str(q))
        fout.write("\t")
        fout.write(str(vals / float(q)))
        fout.write("\n")
    fout.close()

# LOOP THROUGH EACH FILE
for file in files:
    if "_50_moving_25bp" in file and ".xlsx" not in file:
        print (file)
        name, rest = file.split("_50_moving_25bp")
        getVariable(name)

