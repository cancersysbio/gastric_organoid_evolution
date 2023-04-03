#!/usr/bin/env python
# coding: utf-8


## GENERAL ABBREVIATIONS ##
# bc = 10X BarCode
# umi = unique molecular identifier
# adt = antibody-derived tags
# name =  adt-associated samplename


import os
import gzip
import sys
import operator
from umi_tools.network import UMIClusterer
clusterer = UMIClusterer(cluster_method="directional")
import numpy as np



polyT = 'TTTTTTTTTT'
seq_threshold = 27 # for 0.2% error rate



## Function: hamming ##
    # to account for sequencing error or somatic mutation
    # operator.ne return TRUE if str1 != str2
        # 'sum' function return the number of difference
def hamming1(string1, string2):
    '''
    Compare two strings if they have the same length.
    Return number of difference among the two strings.
    '''
    if len(string1) == len(string2):
        return sum(map(operator.ne, string1, string2)) 



## Function: extract 10x barcode (bc), UMI, and cell hashing barcode ##
# out: dictionary with bc_umi as the key and {ecb1:count, ecb2:count, ...} as the value

# FastQ file: each entry has 4 lines
    # Line 1: sequencer info and read ID 
        # e.g. @M02357:509:000000000-CJ9TR:1:1101:19639:1975 1:N:0:ATTACTCG+TCTTTCCC
    # Line 2: read
    # Line 3: +/- strands
    # Line 4: quality check

def extract_fastQ(fq_dir, fastqname, hashhandle, hamming_distance):

    R1_name = fastqname + 'R1.fq.gz' # e.g. 'output_HASH_ECB-ATTACTCG-R1.fq.gz'
    R2_name = fastqname + 'R2.fq.gz' # e.g. 'output_HASH_ECB-ATTACTCG-R2.fq.gz'

    if not os.path.exists(os.path.join(fq_dir, R1_name)) and not os.path.exists(os.path.join(fq_dir, R2_name)):
        R1_name = fastqname + 'R1.fastq.gz' # e.g. 'output_HASH_ECB-ATTACTCG-R1.fastq.gz'
        R2_name = fastqname + 'R2.fastq.gz' # e.g. 'output_HASH_ECB-ATTACTCG-R2.fastq.gz'
    
    if not os.path.exists(os.path.join(fq_dir, R1_name)) and not os.path.exists(os.path.join(fq_dir, R2_name)):
        R1_name = fastqname + 'R1.txt.gz' # e.g. 'output_HASH_ECB-ATTACTCG-R1.txt.gz'
        R2_name = fastqname + 'R2.txt.gz' # e.g. 'output_HASH_ECB-ATTACTCG-R2.txt.gz'

    if not os.path.exists(os.path.join(fq_dir, R1_name)) and not os.path.exists(os.path.join(fq_dir, R2_name)):
        print('Error: input fastq file format must be .fq.gz or .fastq.gz or .txt.gz')

    # R1: line 2 read = 16+10+10 nt 
    # 10x barcode (bc)(16 nt) and UMI (10 nt), followed by the polyT seq (10 nt)
    # R2: line 2 read = 20 nt ebc handle (customized), followed by 30 nt ebc
    temp_dict = {}
    R1R2_dict = {}
    adt_ham_list = []
    
    counter = 0
    with gzip.open(os.path.join(fq_dir, R1_name), 'rt') as R1_file, gzip.open(os.path.join(fq_dir, R2_name), 'rt') as R2_file:
        for line1, line2 in zip(R1_file, R2_file): # parse R1 and R2 fastQ at the same time
            line1 = line1.strip()             
            
            line2 = line2.strip()
            # Extract read ID: read line 1 of each entry
            if counter % 4 == 0: # line 1, line 5, line 9, ...
                readID1 = ':'.join(line1.split(' ')[0].split(':')[4:]) # e.g. 1101:19639:1975
                readID2 = ':'.join(line2.split(' ')[0].split(':')[4:])

            # read every 2nd line: actual read
            if counter % 4 == 1: # line 2, line 6, line 10, ... 
                # To prevent truncated reads from creating out of range error
                if len(line1) >= 36: # 36 is the length of bc, umi, and polyT.
                    polyT_ham = hamming1(line1[26:36], polyT) # test polyT sequence (quality control)
                    if polyT_ham <= hamming_distance: # if quality is good enough
                        bc_umi = line1[:26]
                        if len(line2) >= 50:
                            hash_handle_ham = hamming1(line2[0:21], hashhandle) # test ecb pcr handle quality
                            if hash_handle_ham <= hamming_distance:
                                ADT = line2[21:33]
                                if readID1 == readID2:
                                    if bc_umi not in R1R2_dict:
                                        R1R2_dict[bc_umi] = {ADT:1} # store bc_UMI and ecb
                                    else:
                                        temp_dict = R1R2_dict[bc_umi] 
                                        if ADT not in temp_dict:
                                            temp_dict[ADT] = 1
                                        else:
                                            temp_dict[ADT] += 1 # count how many when more than 1
                                        R1R2_dict[bc_umi] = temp_dict
                                    adt_ham_list.append(ADT)

                                else:
                                    print('%s is not the same as %s!!!'% (readID1, readID2))
            
                        
            counter += 1
            
            if counter%500000 == 0:
                print('Finished %d lines' % counter)
                
    return R1R2_dict, adt_ham_list




## Functions: for bc hamming ##

def convertForw(str16): 
    return ['ACGTN'.index(ch) for ch in str16]

def convertBack(lst):
    return ''.join(['ACGTN'[i] for i in lst])

def make_bc10x(bc_dir, bcfile):
    ### Get all 10X defined cell barcodes
    with open(os.path.join(bc_dir, bcfile), "r") as fd:
        return [convertForw(line.split()[0][:-2]) for line in fd.readlines()]  

def make_adt(adtref):
    ### Get all adt reference
    return [convertForw(a) for a in adtref]   

def hamming2(bc,bc10x): # 'ACTGTGCACACTGTAC', list of list
    bc = convertForw(bc)
    bc = np.array(bc)
    diff = bc - bc10x
    return np.count_nonzero(diff, axis=1)

def countOnes(arr): 
    return np.count_nonzero(arr == 1)

def testHamming(a, b, expected):
    b = [convertForw(item) for item in b]
    result = hamming2(a,b)
    for i in range(len(a)):
        if expected[i] != result[i]: return False
    return True

def handleOne(bc, bc10x):
    arrHamming = hamming2(bc, bc10x)
    mindex = np.argmin(arrHamming)
    minHam = arrHamming[mindex]
    b = convertBack(bc10x[mindex])
    if minHam == 0: return b
    if minHam == 1:
        if countOnes(arrHamming) == 1: return b  # if unique
        return "Remove_multham1"
    return "Remove_Minham"

def handleTwo(bc, bc10x):
    arrHamming = hamming2(bc, bc10x)
    mindex = np.argmin(arrHamming)
    minHam = arrHamming[mindex]
    b = convertBack(bc10x[mindex])
    if minHam < 3: return b
    if minHam > 2:
        if countOnes(arrHamming) == 1: return b  # if unique
        return "Remove_multham1"
    return "Remove_Minham"

def testHandleOne(a,b):
    b = [convertForw(item) for item in b]
    return handleOne(a,b)


def testHandleTwo(a,b):
    b = [convertForw(item) for item in b]
    return handleTwo(a,b)



def count(arr,n): return np.count_nonzero(arr == n)

def Hamming_fast(a,b): # 'ACTGTGCACACTGTAC', list of list
    b = [convertForw(item) for item in b]
    result = hamming2(a,b)
    return (result)

# print (chrisHamming('ACGT',['ACGT','AAGT','AAAT','AAAA','TGCA']))





## Tests for the above functions ##

# assert convertForw('ACGT') == [0,1,2,3]
# assert convertForw('AAAA') == [0,0,0,0]
# assert convertForw('TTTT') == [3,3,3,3]

# assert convertBack([0,1,2,3]) == 'ACGT'
# assert convertBack([0,0,0,0]) == 'AAAA'
# assert convertBack([3,3,3,3]) == 'TTTT'

# assert testHamming('ACGT',['ACGT','AAGT','AAAT','AAAA','TGCA'],[0,1,2,3,4])

# assert testHandleOne('ACGT',['ACGT','AAGT','AAAT','AAAA','TGCA']) == 'ACGT'
# assert testHandleOne('ACGT',['ACGT','ACGT']) == 'ACGT'
# assert testHandleOne('ACGT',['AAGT','AAAT','AAAA','TGCA']) == 'AAGT'
# assert testHandleOne('ACGT',['AAAT','AAAA','TGCA']) == "Remove_Minham"
# assert testHandleOne('ACGT',['AAGT','CCGT','AAAA','TGCA']) == "Remove_multham1"
# assert testHandleOne('ACGT',['AAGT','AAGT']) == "Remove_multham1"






