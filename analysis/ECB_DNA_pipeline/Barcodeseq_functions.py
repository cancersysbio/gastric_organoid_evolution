import os
import gzip

#### THIS FUNCTION EXTRACTS ECB READS FROM FASTQ FILE
def extract_ECB(name,curdir,basedir,ecb_handle_1,ecb_handle_2,fout_sum):
    print (name)
    os.chdir(curdir)
    if not os.path.exists("ECB"):
        os.makedirs("ECB")

    fout = open("ECB/"+name + "_ECB.txt", "w")  # OUTPUT FILES CONTAINING ECBS IN COLUMN 1 AND NUMBER READS PER ECB IN COLUMN 2

    c = 0  # LINE COUNTER
    n = 0  # BARCODE COUNTER
    q = 0  # READS WITH NO ECB COUNTER
    totreads = 0
    readCountDic = {}  # DICTIONARY CONTAINING ECBS AS A KEYS AND NUMBER OF OCCURANCES AS VALUES

    fd = gzip.open(basedir + "/fastqs/" + name + ".fastq.gz", "r")
    line = fd.readline().decode()
    print (line)
    while line:
        line = line.split("\t")
        if c % 4 == 1:                                                                  # LINE WITH SEQUENCES IN FASTQ FILE
            #print(line[0][0:25])
            totreads += 1
            if line[0][0:25] == ecb_handle_1 and line[0][55:76] == ecb_handle_2:        # WERE ONLY KEEPING READS WITH CORRECT HANDLES
                n += 1
                ct = line[0][25:55]                                                     # THIS IS WHERE THE ECB BARCODE IS LOCATED
                if "N" not in ct:
                    if ct in readCountDic:
                        readCountDic[ct] += 1
                    if ct not in readCountDic:
                        readCountDic[ct] = 1
            if line[0][0:25] != ecb_handle_1 or line[0][55:76] != ecb_handle_2:
                q += 1                                                                  # IF HANDLES DONT HAVE THE CORRECT SEQUENCE, WE DISCARD THE READ
        c += 1
        line = fd.readline().decode()
    fout_sum.write(str(name))
    fout_sum.write("\t")
    fout_sum.write(str(totreads))
    fout_sum.write("\t")
    fout_sum.write(str(n))
    fout_sum.write("\t")
    fout_sum.write(str(q))
    fout_sum.write("\n")

    for d in readCountDic:
        fout.write(str(d[0:30]))
        fout.write("\t")
        fout.write(str(readCountDic[d]))
        fout.write("\n")


#### THIS FUNCTION ADDS ECBs FROM NEW FILE TO A GROWING DICTIONARY OF ECBs
def ecbdic(fd,dic):  # INPUT NEW FILE WITH ECBS AND COUNTS PER ECB, AND A DICTIONARY CONTAINING ECBS AND READ COUNTS PER ECB
    line = fd.readline()
    while line:
        line = line.split()
        ecb = line[0]
        reads = int(line[1])
        if ecb in dic:            # JUST ADD ECB AND READS IF ECB IN DIC
            cur_reads = dic[ecb]
            new_reads = cur_reads + int(reads)
            dic[ecb] = new_reads
        if ecb not in dic:        # ADD ECB AND READS IF ECB NOT IN DIC
            dic[ecb] = reads
        line = fd.readline()
    return dic


#### THIS FUNCTION PRINTS OUT READ GROUPS PER SAMPLE

def RG_loop(out_name,master,RG_folder):
    if not os.path.exists(RG_folder):
        os.makedirs(RG_folder)

    fdp = open("ECB/"+out_name + "_ECB.txt", "r")
    fout = open(RG_folder+"/"+out_name + "_RG.txt", "w")

    line = fdp.readline()
    dic = {}
    while line:
        line = line.split()
        ecb = line[0]
        reads = int(line[1])
        if ecb in master:
            readgroup = master[ecb]
            if readgroup in dic:
                current_reads = dic[readgroup]
                upd_reads = current_reads+reads
                dic[readgroup] = upd_reads
            if readgroup not in dic:
                dic[readgroup] = reads

        if ecb not in master:
            if "Na" in dic:
                current_reads = dic["Na"]
                upd_reads = current_reads+reads
                dic["Na"] = upd_reads
            if "Na" not in dic:
                dic["Na"] = reads

        line = fdp.readline()

    for d in dic:
        fout.write(str(d))
        fout.write("\t")
        fout.write(str(dic[d]))
        fout.write("\n")

### FUNCTION TO CREATE READ GROUP MATRIX

def rgcMat(names, fout, RG_folder):
    RGS = []

    # PUT RGS FROM ALL SAMPLES INTO THE LIST RGS
    for n in names:
        fd = open(RG_folder+"/" + n + "_RG.txt", "rU")
        line = fd.readline()
        while line:
            line = line.split()
            RG = line[0]
            if RG != "Na":
                RG = int(RG)
                if RG not in RGS:
                    RGS.append(RG)
            line = fd.readline()
        fd.close()

    RGS.sort()

    fout.write("sample")
    for r in RGS:
        fout.write("\t")
        fout.write("RG"+str(r))
    fout.write("\n")

    for n in names:
        rgdic = {}
        fd = open(RG_folder+"/" + n + "_RG.txt", "rU")
        line = fd.readline()
        while line:
            line = line.split()
            RG = line[0]
            count = int(line[1])
            rgdic[RG] = count
            line = fd.readline()
        fout.write(str(n))
        for r in RGS:
            r = str(r)
            fout.write("\t")
            if r in rgdic:
                fout.write(str(rgdic[r]))
            if r not in rgdic:
                fout.write("0")
        fout.write("\n")
        fd.close()



### FUNCTIONS BELOW ARE TO CREATE MULLER PLOTS

def get_totsum(fd):
    line = fd.readline()
    totsum = 0
    while line:
        line = line.split()
        if line[0] != "Na":
            totsum += int(line[1])
        line = fd.readline()
    fd.close()
    return totsum

def get_all_RG(fd,all_RGs,totsum):
    line = fd.readline()
    while line:
        line = line.split()
        readgroup = line[0]
        count = float(line[1])
        if count / float(totsum+0.0000000001) >= 0.001:
            if readgroup != "Na":
                all_RGs[readgroup] = 0
        line = fd.readline()
    return all_RGs

def get_RG(out_name,RG_folder):
    fd = open(RG_folder+"/"+out_name+"_RG.txt", "rU")
    RGs = {}
    newRGs = {}
    line = fd.readline()
    totsum = 0
    while line:
        line = line.split()
        readgroup = line[0]
        count = int(line[1])
        if readgroup != "Na":
            RGs[readgroup] = count
            totsum += count
        line = fd.readline()
    fd.close()
    fd = open(RG_folder+"/"+out_name+"_RG.txt", "rU")
    line = fd.readline()
    newTotsum = 0
    while line:
        line = line.split()
        readgroup = line[0]
        count = int(line[1])
        if readgroup != "Na":
            if count / float(totsum+0.0000000001) >= 0.001:
                newTotsum += count
                newRGs[readgroup] = count
        line = fd.readline()
    fd.close()
    return RGs,newRGs,totsum,newTotsum

def printMuller(fout,code,out_names,RG_folder):
    fout.write("Generation")
    fout.write("\t")
    fout.write("Identity")
    fout.write("\t")
    fout.write("Population")
    fout.write("\t")
    fout.write("Frequency")
    fout.write("\t")
    fout.write("Group_id")
    fout.write("\t")
    fout.write("Unique_id")
    fout.write("\t")
    fout.write("RG")
    fout.write("\t")
    fout.write("TimePoint")
    fout.write("\n")
    generation = 1
    all_RGs = {}
    for out_name in out_names:
        if "parent" in out_name or "Parent" in out_name or code in out_name:
            fd = open(RG_folder+"/"+out_name + "_RG.txt", "rU")
            totsum = get_totsum(fd)
            fd.close()
            fd = open(RG_folder+"/"+out_name + "_RG.txt", "rU")
            all_RGs = get_all_RG(fd, all_RGs, totsum)
            fd.close()
    for out_name in out_names:
        if "parent" in out_name or "Parent" in out_name or code in out_name:
            RGs,newRGs,totsum,newTotsum = get_RG(out_name,RG_folder)
            for r in all_RGs:
                curclone = "clone_" + str(r)
                if r in newRGs:
                    fout.write(str(generation))
                    fout.write("\t")
                    fout.write(curclone)
                    fout.write("\t")
                    fout.write(str(RGs[r]))
                    fout.write("\t")
                    fout.write(str(float(RGs[r])/float(newTotsum)))
                    fout.write("\t")
                    fout.write(curclone)
                    fout.write("\t")
                    fout.write(curclone+"_"+str(generation))
                    fout.write("\t")
                    fout.write(str(r))
                    fout.write("\t")
                    fout.write(str(out_name))
                    fout.write("\n")
                if r not in newRGs:
                    fout.write(str(generation))
                    fout.write("\t")
                    fout.write(curclone)
                    fout.write("\t")
                    fout.write("0")
                    fout.write("\t")
                    fout.write("0")
                    fout.write("\t")
                    fout.write(curclone)
                    fout.write("\t")
                    fout.write(curclone+"_"+str(generation))
                    fout.write("\t")
                    fout.write(str(r))
                    fout.write("\t")
                    fout.write(str(out_name))
                    fout.write("\n")
            generation+=1

