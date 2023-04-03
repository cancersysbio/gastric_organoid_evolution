#############################################################################################################################
##
##  SUMMARY FRACTION GENOME ALTERED
##
##  SUMMARIZES FRACTION GENOME ALTERED PER CULTURE, USED IN FIGURE 1D
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################

from statistics import mean
import os

### SET PATH
path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_1/"
os.chdir(path)

### GET ALL FILE NAMES
files = os.listdir(path+"INPUT/CNA_smooth_50kb_shallow_MA25_Centered_1/Summary")

### CHECK IF PATH FOR NEW FILES EXISTS, OTHERWISE CREATE

# Check whether the specified path exists or not
isExist = os.path.exists("PLOTS")
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs("PLOTS")

isExist = os.path.exists("PLOTS/Figure_1D_Fraction_genome_altered")
if not isExist:
    os.makedirs("PLOTS/Figure_1D_Fraction_genome_altered")

### OPEN FILES TO WRITE IN
# FRACTION GENOME ALTERED, AVERAGE PER WINDOW
fout = open("PLOTS/Figure_1D_Fraction_genome_altered/Summary_fraction_genome_altered_per_win_movingAve_25p.csv", "w")
# FRACTION GENOME ALTERED, AVERAGE ACROSS CHROMOSOME (NORMALIZED FOR CHROMOSOME LENGTH
fout2 = open("PLOTS/Figure_1D_Fraction_genome_altered/Summary_fraction_genome_altered_per_chr_movingAve_25p.csv", "w")

fout.write("Sample"+","+"Donor"+","+"Clone"+","+"ECB_sample"+","+"Day"+","+"total_nr_win"+","+"total_nr_win_up"+","+"total_nr_win_down"+","+
           "fraction_win_up"+","+"fraction_win_down"+","+"fraction_total_win_altered"+"\n")
fout2.write("Sample"+","+"Donor"+","+"Clone"+","+"ECB_sample"+","+"Day"+","+"average_fraction_win_up"+","+"average_fraction_win_down"+","+
            "average_fraction_total_win_altered"+"\n")

def getAve(file,name,fout,fout2):
    fd = open("INPUT/CNA_smooth_50kb_shallow_MA25_Centered_1/Summary/"+file, "r")
    line = fd.readline()
    line = fd.readline()

    ### PARAMETERS
    total_nr_win = 0
    total_nr_win_up = 0
    total_nr_win_down = 0
    ecb = "NO"

    chrarm_nr_percent_up = []
    chrarm_nr_percent_down = []
    chrarm_nr_percent_total_altered = []

    if "D1" in name:
        donor = "D1"
    if "D2" in name:
        donor = "D2"
    if "D3" in name:
        donor = "D3"
    if "WT" in name:
        clone = "WT"
    if "C1" in name:
        clone = "C1"
    if "C2" in name:
        clone = "C2"
    if "C3" in name:
        clone = "C3"

    day = name[5:-1]

    ### ADJUST DAYS DEPENDING ON SAMPLE
    if "WT" in name:
        day = name[5:]
    if "Parent" in name:
        day = "Parent"
        ecb = "YES"
    if "R" in name:
        day = name[7:-1]
        ecb = "YES"

    fout.write(str(name)+","+donor+","+clone+","+ecb+","+str(day)+",")
    fout2.write(str(name)+","+donor+","+clone+","+ecb+","+str(day)+",")


    while line:
        line = line.split()
        chrarm = line[0]
        if chrarm != "chr21p": ### CHR21 has 0 windows after moving average
            print (line)
            print (line[3])

            nr_win = int(line[3])
            nr_win_up = int(line[4])
            nr_win_down = int(line[5])

            total_nr_win += nr_win
            total_nr_win_up += nr_win_up
            total_nr_win_down += nr_win_down

            ### CALCULATE PERCENT OF WINDOWS THAT ARE ALTERED
            chrarm_nr_percent_up.append(float(nr_win_up)/float(nr_win))
            chrarm_nr_percent_down.append(float(nr_win_down)/float(nr_win))
            chrarm_nr_percent_total_altered.append(float(nr_win_up+nr_win_down)/float(nr_win))
        line = fd.readline()

    ### WRITE TO FILE
    fout.write(str(total_nr_win)+","+str(total_nr_win_up)+","+str(total_nr_win_down)+","+
               str(float(total_nr_win_up)/float(total_nr_win))+","+str(float(total_nr_win_down)/float(total_nr_win))+","+
               str(float(total_nr_win_up+total_nr_win_down)/float(total_nr_win))+"\n")
    fout2.write(str(name)+","+str(mean(chrarm_nr_percent_up))+","+str(mean(chrarm_nr_percent_down))+","+
                str(mean(chrarm_nr_percent_total_altered))+"\n")
    fd.close()

### LOOP THROUGH ALL FILES
for file in files:
    if "_moving_25bp_average" in file:
        name,rest = file.split("_50_moving_25bp_average")
        getAve(file,name,fout,fout2)

