#############################################################################################################################
##
##  EXTRACT MUTATIONAL INFORMATION FROM .MAF FILES
##  Get nr SNV, Coding SVN, SNV in Cosmic, Indels, Location
##
##  (for supplementary table 5 and figures 2A and S5B)
##
##
##  Author: Kasper Karlsson
##
##
############################################################################################################################


### Get nr SNV, Coding SVN, SNV in Cosmic, Indels, Location (for supplementary table 5 and figures 2A and S5B)

import os

path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/Mutect2_multiallele"
os.chdir(path)

fout1 = open("Summary_table_Mutations.txt","w")
fout2 = open("Summary_table_codingAlt.txt","w")
fout3 = open("Summary_table_codingAlt_cosmic.txt","w")
fout4 = open("Summary_table_indels.txt","w")
fout5 = open("Summary_table_location.txt","w")
fout6 = open("Summary_table_cosmic.txt","w")


fd = open("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/sequtils/Supporting_files/cancer_gene_census_Cosmic_v94.csv","r")
line = fd.readline()
cosmic_muts = {}
while line:
    line = line.split(",")
    gene = line[0]
    if gene not in cosmic_muts:
        cosmic_muts[gene] = 0
    line = fd.readline()

print(cosmic_muts)

fd = open("All_mutations_combined.txt","r")
line = fd.readline()

snvDic = {}
codingAlts = {}
codingAlts_cosmic = {}
nrDel = {}
nrIns = {}
locDic = {}
cosmic = {}


while line:
    line = line.split()
    sample = line[0]
    gene = line[1]
    location = line[9]
    type = line[10]

    if location != "Variant_Classification":
        coding_locations = ["Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins"]

        if type == "SNP":
            if sample in snvDic:
                snvDic[sample] += 1
            if sample not in snvDic:
                snvDic[sample] = 1
        if location in coding_locations:
            if sample in codingAlts:
                codingAlts[sample] += 1
            if sample not in codingAlts:
                codingAlts[sample] = 1
                if gene in cosmic_muts:
                    if sample in codingAlts_cosmic:
                        codingAlts_cosmic[sample] += 1
                    if sample not in codingAlts_cosmic:
                        codingAlts_cosmic[sample] = 1
        if type == "DEL":
            if sample in nrDel:
                nrDel[sample] += 1
            if sample not in nrDel:
                nrDel[sample] = 1
        if type == "INS":
            if sample in nrIns:
                nrIns[sample] += 1
            if sample not in nrIns:
                nrIns[sample] = 1


        if location == "Missense_Mutation" or location == "Nonsense_Mutation" or location == "In_Frame_Del" or location == "Frame_Shift_Ins" or location == "Frame_Shift_Del" or location == "Silent" or location == "In_Frame_Ins":
            location = "Coding"
        if location == "Splice_Region" or location == "Splice_Site":
            location = "Splice_Site"
        if location == "3'Flank" or location == "5'Flank" or location == "IGR":
            location = "Intergenic"


        sampleLoc = sample+":"+location


        if sampleLoc in locDic:
            locDic[sampleLoc] += 1
        if sampleLoc not in locDic:
            locDic[sampleLoc] = 1

        if gene in cosmic_muts:
            if sample in cosmic:
                cosmic[sample] += 1
            if sample not in cosmic:
                cosmic[sample] = 1


    line = fd.readline()


print (snvDic)

fout1.write("Sample"+"\t"+"Patient"+"\t"+"Culture"+"\t"+"Annotation"+"\t"+"Timepoint"+"\t"+"Count"+"\n")
for sample in snvDic:
    if "MID" in sample or "LATE" in sample:
        DC,days,timepoint = sample.split("_")
        anno = "EML"
    if "MID" not in sample and "LATE" not in sample:
        DC,timepoint = sample.split("_")
        anno = "Timecourse"
    if "D2C2_729d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "729d"
    if "D3C1_718d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "718d"
    donor = DC[0:2]
    culture = DC[2:4]
    nrSNV = snvDic[sample]
    fout1.write(str(sample) + "\t" + str(donor) + "\t" + str(culture) + "\t" + str(anno) + "\t" + str(timepoint) + "\t")
    fout1.write(str(nrSNV) + "\n")

fout2.write("Sample"+"\t"+"Patient"+"\t"+"Culture"+"\t"+"Annotation"+"\t"+"Timepoint"+"\t"+"Count"+"\n")
for sample in codingAlts:
    if "MID" in sample or "LATE" in sample:
        DC,days,timepoint = sample.split("_")
        anno = "EML"
    if "MID" not in sample and "LATE" not in sample:
        DC,timepoint = sample.split("_")
        anno = "Timecourse"
    if "D2C2_729d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "729d"
    if "D3C1_718d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "718d"
    donor = DC[0:2]
    culture = DC[2:4]
    nrCodingAlts = codingAlts[sample]
    fout2.write(str(sample) + "\t" + str(donor) + "\t" + str(culture) + "\t" + str(anno) + "\t" + str(timepoint) + "\t")
    fout2.write(str(nrCodingAlts) + "\n")

fout3.write("Sample"+"\t"+"Patient"+"\t"+"Culture"+"\t"+"Annotation"+"\t"+"Timepoint"+"\t"+"Count"+"\n")

for sample in codingAlts_cosmic:
    if "MID" in sample or "LATE" in sample:
        DC,days,timepoint = sample.split("_")
        anno = "EML"
    if "MID" not in sample and "LATE" not in sample:
        DC,timepoint = sample.split("_")
        anno = "Timecourse"
    if "D2C2_729d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "729d"
    if "D3C1_718d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "718d"
    donor = DC[0:2]
    culture = DC[2:4]
    nrCodingAlts_cosmic = codingAlts_cosmic[sample]
    fout3.write(str(sample) + "\t" + str(donor) + "\t" + str(culture) + "\t" + str(anno) + "\t" + str(timepoint) + "\t")
    fout3.write(str(nrCodingAlts_cosmic) + "\n")

fout4.write("Sample"+"\t"+"Patient"+"\t"+"Culture"+"\t"+"Annotation"+"\t"+"Timepoint"+"\t"+"Type"+"\t"+"Count"+"\n")
for sample in nrDel:
    if "MID" in sample or "LATE" in sample:
        DC,days,timepoint = sample.split("_")
        anno = "EML"
    if "MID" not in sample and "LATE" not in sample:
        DC,timepoint = sample.split("_")
        anno = "Timecourse"
    if "D2C2_729d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "729d"
    if "D3C1_718d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "718d"
    donor = DC[0:2]
    culture = DC[2:4]
    nrDels = nrDel[sample]
    type = "DEL"
    fout4.write(str(sample) + "\t" + str(donor) + "\t" + str(culture) + "\t" + str(anno) + "\t" + str(timepoint) + "\t")
    fout4.write(str(type) + "\t" + str(nrDels) + "\n")

for sample in nrIns:
    if "MID" in sample or "LATE" in sample:
        DC,days,timepoint = sample.split("_")
        anno = "EML"
    if "MID" not in sample and "LATE" not in sample:
        DC,timepoint = sample.split("_")
        anno = "Timecourse"
    if "D2C2_729d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "729d"
    if "D3C1_718d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "718d"
    donor = DC[0:2]
    culture = DC[2:4]
    nrInss = nrIns[sample]
    type = "INS"
    fout4.write(str(sample) + "\t" + str(donor) + "\t" + str(culture) + "\t" + str(anno) + "\t" + str(timepoint) + "\t")
    fout4.write(str(type) + "\t" + str(nrInss) + "\n")

sumLocPerSample = {}
for sampleLoc in locDic:
    sample, location = sampleLoc.split(":")
    nrAlts = locDic[sampleLoc]
    if sample in sumLocPerSample:
        sumLocPerSample[sample] += nrAlts
    if sample not in sumLocPerSample:
        sumLocPerSample[sample] = nrAlts


fout5.write("Sample"+"\t"+"Patient"+"\t"+"Culture"+"\t"+"Annotation"+"\t"+"Timepoint"+"\t"+"Type"+"\t"+"Count"+"\t"+"Fraction"+"\n")
for sampleLoc in locDic:
    sample,location = sampleLoc.split(":")
    if "MID" in sample or "LATE" in sample:
        DC,days,timepoint = sample.split("_")
        anno = "EML"
    if "MID" not in sample and "LATE" not in sample:
        DC,timepoint = sample.split("_")
        anno = "Timecourse"
    if "D2C2_729d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "729d"
    if "D3C1_718d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "718d"
    donor = DC[0:2]
    culture = DC[2:4]
    count = locDic[sampleLoc]
    sumAlts = sumLocPerSample[sample]
    fraction = float(count) / float(sumAlts)
    type = location
    fout5.write(str(sample) + "\t" + str(donor) + "\t" + str(culture) + "\t" + str(anno) + "\t" + str(timepoint) + "\t")
    fout5.write(str(type) + "\t" + str(count) + "\t" + str(fraction) +"\n")


fout6.write("Sample"+"\t"+"Patient"+"\t"+"Culture"+"\t"+"Annotation"+"\t"+"Timepoint"+"\t"+"Count"+"\n")
for sample in cosmic:
    if "MID" in sample or "LATE" in sample:
        DC,days,timepoint = sample.split("_")
        anno = "EML"
    if "MID" not in sample and "LATE" not in sample:
        DC,timepoint = sample.split("_")
        anno = "Timecourse"
    if "D2C2_729d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "729d"
    if "D3C1_718d_LATE" in sample:
        anno = "Timecourse"
        timepoint = "718d"
    donor = DC[0:2]
    culture = DC[2:4]
    nrCosmic = cosmic[sample]
    fout6.write(str(sample) + "\t" + str(donor) + "\t" + str(culture) + "\t" + str(anno) + "\t" + str(timepoint) + "\t")
    fout6.write(str(nrCosmic) + "\n")