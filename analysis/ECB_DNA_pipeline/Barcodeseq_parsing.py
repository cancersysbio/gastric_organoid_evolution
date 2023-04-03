import os
import sys
from barcodeseq_functions import *
from umi_tools.network import UMIClusterer
clusterer = UMIClusterer(cluster_method = 'directional')

### Specify input file

basedir = sys.argv[1]   # DIRECTORY WITH ALL SAMPLES TO RUN ECB EXTRACTION ON
subdir = sys.argv[2]    # SUB DIRECTORY THAT SHOULD CONTAIN SAMPLE NAMEFILE (E.G. "T12" TO RUN A SUBSET OF TIME POINTS)
namefile = sys.argv[3]  # INPUT FILE WIHT SAMPLENAMES WITHOUT fastq.gz
masterfile = sys.argv[4]  # INPUT ECB DECONVOLVED FILE (STANDARD = ECB_groups_master_list.txt)
RG_folder = sys.argv[5] # READ GROUP FOLDER (USE E.G. "BASE" OR "PAIRS" IF MASTERLIST HAS BEEN UPDATED WITH DUAL INSERTION EVENTS)

##----------------------Variables------------------------------##
ecb_handle_1 = "ACTGACTGCAGTCTGAGTCTGACAG"
ecb_handle_2 = "AGCAGAGCTACGCACTCTATG"


### READ IN ALL FILE NAMES INTO A LIST
curdir = basedir+"/"+subdir
print (curdir)
os.chdir(curdir)

print ("namefile",basedir+"/"+subdir+"/"+namefile)
fd = open(namefile+".txt","r")
line = fd.readline()
names = []
while line:
    line = line.split()
    names.append(line[0])
    line = fd.readline()

### EXTRACT AND OUTPUT ECBs FROM EACH FILE
print ("Extracting ECBs")

fout_sum = open("Summary_barcodes_per_file.txt","w")
fout_sum.write("sample\ttotal_reads\tnumber_ok_barcodes\tnumber_wrong_reads\n")

for name in names:
    extract_ECB(name,curdir,basedir,ecb_handle_1,ecb_handle_2,fout_sum)

fout_sum.close()

### COMBINE ALL ECBs TO ONE FILE
print ("Combining ECBs")

fout = open("ECB_all_files.txt","w")

comb = {}
for name in names:
    fd = open("ECB/"+name + "_ECB.txt")
    comb = ecbdic(fd,comb)
    fd.close()

for c in comb:
    if int(comb[c]) >= 2:     # ONLY PRINT OUT ECBS WITH MORE THAN ONE COUNT, REDUCES NOICE
        fout.write(c)
        fout.write("\t")
        fout.write(str(comb[c]))
        fout.write("\n")

fout.close()

### RUN UMITOOLS TO ASSOCIATE BARCODES WITH READ GROUPS

print ("Running UMItools Clusterer")

fout = open(RG_folder+"/ECB_groups_master_list.txt", "w")
fd = open("ECB_all_files.txt","r")
line = fd.readline()
umidic={}
while line:                         # Create a dictionary with ECBs and corresponding read counts as input for the clusterer
    line=line.split()
    umi=line[0]
    umi_as_bytes=str.encode(umi)
    umidic[umi_as_bytes]=int(line[1])
    line = fd.readline()

final_umis = clusterer(umidic, 1)
n=0                                 # RG number, also counts total number of read groups
totBarcodes = 0                     # Count total number of different barcodes

for f in final_umis:
    for each in f:
        fout.write(str(each.decode()))
        fout.write("\t")
        fout.write(str(n))
        fout.write("\n")
        totBarcodes += 1
    n+=1

fout.close()
fout2 = open("Summary_Clustering.txt", "w")
fout2.write("totBarcodes")
fout2.write("\t")
fout2.write(str(totBarcodes))
fout2.write("\n")
fout2.write("totReadGroups")
fout2.write("\t")
fout2.write(str(n))
fout2.close()



### READ IN MASTERFILE WITH BARCODES AND ASSOCIATED READ GROUP

fdm = open(masterfile, "r")

line = fdm.readline()
master = {}
while line:
    line = line.split()
    ecb = line[0]
    group = int(line[1])
    master[ecb] = group
    line = fdm.readline()
fdm.close()

### COUNT READS PER GROUP FOR YOUR SAMPLE

for out_name in names:
    RG_loop(out_name,master,RG_folder)


### CREATE RGC MATRIX

fout = open(RG_folder+"/RGC_matrix.txt", "w")

rgcMat(names,fout,RG_folder)
fout.close()



### CREATE FILES FOR MULLER PLOTS

fout1 = open(RG_folder+"/Mueller_R1_df.txt","w")
fout2 = open(RG_folder+"/Mueller_R2_df.txt","w")
fout3 = open(RG_folder+"/Mueller_R3_df.txt","w")

printMuller(fout1,"R1",names,RG_folder)
printMuller(fout2,"R2",names,RG_folder)
printMuller(fout3,"R3",names,RG_folder)

print ("done")
