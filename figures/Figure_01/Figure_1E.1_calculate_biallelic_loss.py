#############################################################################################################################
##
##  CALCULATE SEGMENTS OF BI-ALLELIC LOSS
##
##  SEGMENTS ARE CALLED AS LOST IF A GENOMIC WINDOW CONTAINS 25% OR LESS READS COMPARED TO MEDIAN NUMBER OF READS PER WINDOW
##
##  OUTPUT IS USED FOR FIGURE 1E
##  FHIT hg38: chr3:59,747,277-61,251,459
##  CDKN2A hg38: chr9:21,967,752-21,995,324
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

import os

### CHECK IF PATH FOR NEW FILES EXISTS, OTHERWISE CREATE
path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_1/"
os.chdir(path)
# Check whether the specified path exists or not
isExist = os.path.exists("PLOTS/Figure_1E_Timing_CNA")
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs("PLOTS/Figure_1E_Timing_CNA")

isExist = os.path.exists("INPUT/CNA_smooth_50kb_shallow/Dual_KO_0.25")
if not isExist:
    os.makedirs("INPUT/CNA_smooth_50kb_shallow/Dual_KO_0.25")


### PRINT SUMMARY OF FHIT AND CDKN2A BIALLELIC LOSS TO FILE
fout2 = open("PLOTS/Figure_1E_Timing_CNA/DualKO_summary_50kb_0.25.csv","w")
fout2.write("Sample" + "," + "Donor" + "," + "Clone" + "," + "Day" + "," + "Aberration" + "\n")

### GET ALL FILE NAMES
path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_1/INPUT/CNA_smooth_50kb_shallow"
os.chdir(path)
files = os.listdir(path)


### FUNCTION TO EXTRACT SEGMENTS WITH BIALLELIC LOSS
FhitDic = {}
Cdkn2aDic = {}

def getDualKO(file,name,FhitDic,Cdkn2aDic):
    fout = open("Dual_KO_0.25/"+name+"_DualKO.txt", "w") # STORE ALL LOST GENOMIC SEGMENTS PER CULTURE AND TIME (NOT JUST FOR FHIT AND CDKN2A)
    fd = open(file, "r")
    line = fd.readline()
    line = fd.readline()

    while line:
        line = line.split()
        value = float(line[4]) # LOG 2 VALUES FROM QDNASEQ OUTPUT
        centerVal = 2**value # CONVERT TO COPY NUMBER
        if  centerVal <= 0.25:  # KEEP ONLY WINDOWS WITH COPY NUMBER LESS THAN 0.25 (I.E. ONE ALLELE LOST AND THE OTHER ALLELE LOST IN 50% OR MORE OF THE POPULATION)
            for l in line:
                fout.write(str(l)+"\t")
            fout.write(str(centerVal)+"\t")

            chr = line[1]
            start = int(line[2])
            end = int(line[3])

            # STORE START POSITION OF WINDOWS IN THE FHIT GENE
            goFHIT = "no"
            if chr == "3":
                if start >= 59747277 and start <= 61251459:
                    goFHIT = "OK"
                if end >= 59747277 and end <= 61251459:
                    goFHIT = "OK"
            if goFHIT == "OK":
                fout.write("FHIT")
                if name not in FhitDic:
                    FhitDic[name] = start

            # STORE START POSITION OF WINDOWS IN THE CDKN2A GENE
            goCDKN2A = "no"
            if chr == "9":
                if start >= 21950000 and start <= 22000001:
                    goCDKN2A = "OK"
                if end >= 21950000 and end <= 22000001:
                    goCDKN2A = "OK"
            if goCDKN2A == "OK":
                fout.write("CDKN2A")
                if name not in Cdkn2aDic:
                    Cdkn2aDic[name] = start
            if goFHIT == "no" and goCDKN2A == "no":
                fout.write("Na")
            fout.write("\n")
        line = fd.readline()
    fd.close()
    fout.close()
    return(FhitDic,Cdkn2aDic)

### LOOP THROUGH ALL FILES
for file in files:
    if "_50_" in file and "_smo_chrarm" in file:
        if "Mid" not in file:
            print (file)
            name, rest = file.split("_50_")
            FhitDic,Cdkn2aDic =getDualKO(file,name,FhitDic,Cdkn2aDic)

### PRINT LOST GENOMIC SEGMENTS FOR FHIT AND CDKN2A TO FILE
for hit in FhitDic:
    donor = hit[0:2]
    clone = hit[2:4]
    day = hit[5:-1]
    if "WT" not in hit and "Parent" not in hit and "R" not in hit:  # Biallelic loss is not investigated for ECB cultures in figure 1
        fout2.write(hit+","+donor+","+clone+","+day+","+"FHIT"+"\n")
for hit in Cdkn2aDic:
    donor = hit[0:2]
    clone = hit[2:4]
    day = hit[5:-1]
    if "WT" not in hit and "Parent" not in hit and "R" not in hit:  # Biallelic loss is not investigated for ECB cultures in figure 1
        fout2.write(hit+","+donor+","+clone+","+day+","+"CDKN2A"+"\n")


