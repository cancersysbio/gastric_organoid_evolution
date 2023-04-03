#!/usr/bin/env python
# coding: utf-8



## GENERAL ABBREVIATIONS ##
# bc = 10X BarCode
# umi = unique molecular identifier
# adt = antibody-derived tags
# name =  adt-associated samplename



import os
import sys
import gzip
from umi_tools.network import UMIClusterer
clusterer = UMIClusterer(cluster_method="directional")
from hash_functions import *
import time
import copy
import operator
import numpy as np
import csv
from collections import Counter
from operator import itemgetter
from datetime import date


##----------------------Variables------------------------------##
input_dir = sys.argv[1]
fastq_name = sys.argv[2]

## Set Paths ##

# Below paths are auto #
output_dir = input_dir + 'hash_outs/'
if not os.path.exists(output_dir):
    os.mkdir(output_dir) 
interpath = output_dir + 'intermediate_output/'
if not os.path.exists(interpath):
    os.mkdir(interpath)

## for fastq parsing ##
hash_handle = 'CCTTGGCACCCGAGAATTCCA'

## reference ##
bc_filename = 'barcodes.tsv'
adt_ref = 'Hash_lookup.txt'

## Some thresholds and ratio for filtering ##
hamming_distance = 2
majority_percent = 50
true_name_threshold = 0.3
rg_cleanup_thres = 2 # used in step 4
umilim = 2 # used in step 4
name_totalcount_lim = 3# used in step 4
major_name_percent = 0.7 # used in step 4, default 0.7

#### Intermediate Output ####  
adt_hammed_outfile = 'Hash_ADT_clusters_expected.txt'
bc_adt_count_outfile = 'Hash_correlation.txt' #Hash_correlation_condensation_clean_min20all_more70p

#### Output ####
seurat_outfile = 'Hash_forSeurat.csv'

##------------------------------------------------------------##


BeginningOfTime = time.time() # Track time for the whole run

## Parse adt fastq ##
if ',' not in fastq_name: # single run
    print('Parsing hash fastq...')
    concat_dict, adt_forhamming_list = extract_fastQ(input_dir, fastq_name, hash_handle, hamming_distance)
    #{bc_umi:{ecb1:count, ecb2:count ...}}
    print('Parsing hash done.')

else: # multiple runs
    c = 0
    fastq_name = fastq_name.split(',')
    for fastq in fastq_name:
        print('Parsing hash %s ...' %fastq)
        c += 1
        temp_concat_dict, temp_adt_forhamming_list = extract_fastQ(input_dir, fastq, hash_handle, hamming_distance)
        print(len(temp_concat_dict))
        if c == 1: # first pair of R1R2
            concat_dict = {**temp_concat_dict}
            # concat_dict = copy.deepcopy(temp_dict) 
            adt_forhamming_list = temp_adt_forhamming_list.copy()
        elif c > 1: # second pair and so on
            for key in temp_concat_dict:
                if key not in concat_dict:
                    concat_dict[key] = temp_concat_dict[key]
                else:
                    subdict = concat_dict[key] # {ecb:read}
                    subdict = Counter(subdict) + Counter(temp_concat_dict[key])
                    concat_dict[key] = subdict
            adt_forhamming_list.extend(temp_adt_forhamming_list)
        print('%s done' %fastq)



## 1. ADT hamming ##
print('adt hamming ...')

# Parse adt ref
hashtags = []
with open(os.path.join(input_dir, adt_ref), 'rt') as file:
    names = {}
    for line in file:
        line = line.split()
        #print (line)
        if line[1] == 'ADT':
            tag_name = line[2]
            adt = line[3].split(",")
            for a in adt:
                names[a] = tag_name
                hashtags.append(a)

#print (hashtags)
# print (names)
# print (names.keys())

adt_all = {}
adt_counts = {}

min_edit_distance = 2

total_count = 0
ambiguous_count = 0
no_match_count = 0

c = 0

for curr_adt in adt_forhamming_list:
    
    if curr_adt in adt_counts:
        adt_counts[curr_adt] += 1
        if adt_all[curr_adt] == "two_matches":
            ambiguous_count += 1
        if adt_all[curr_adt] == "no_match":
            no_match_count += 1

    else:
        adt_counts[curr_adt] = 1

        adt_ham_dist = Hamming_fast(curr_adt, hashtags)
        scores = []
        for each in range(len(hashtags)):
            scores.append([hashtags[each], adt_ham_dist[each]])

        scores = sorted(scores, key=itemgetter(1))

        best_adt, next_adt = scores[0], scores[1]
	#print ("TEST",best_adt, next_adt)
	
        if best_adt[1] == next_adt[1]:
            adt_all[curr_adt] = "two_matches"
            ambiguous_count += 1 # do not add to final dict, cannot determine which adt is better
        elif best_adt[1] > min_edit_distance:
            adt_all[curr_adt] = "no_match"
            no_match_count += 1  # current read adt does not match any in expected list
        else:
            sample_adt = best_adt[0]
            adt_all[curr_adt] = sample_adt
            
    c += 1
    if c%1000000 == 0:
        print('Finished %d reads' %c)


## Intermediate Output ##
# Hash_ADT_clusters_expected.txt

header = ['adt_read', 'adt_counts', 'adt_ref', 'sample_origin']

with open(os.path.join(interpath, adt_hammed_outfile), "w") as output:
    print('\t'.join(header), file=output)
    for adt in adt_all.keys():
        output.write(adt+'\t')
        output.write(str(adt_counts[adt])+'\t')
        output.write(adt_all[adt]+'\t')  # sample adt sequence
        if adt_all[adt] in names:
            output.write(names[adt_all[adt]]+'\n')  # sample name
        else:
            output.write(adt_all[adt]+'\n')

print ("totReads", total_count)
print ("totUniqueADTs",len(adt_all))  #number of ADTs from sequencing data
print ("totClusters",len(names)) #number of clusters output
print ("AmbiguousSampleAssignment", ambiguous_count)
print ("NoSampleAssignment", no_match_count)


# Replace adt with hammed adt. Out = adt_corrected_dict, adt_count_dict
# create a dictionary that pull all bc:
    # bc_umi_adt_dict:   bc:[[umi1, adt1, read1], [umi2, adt2, read2]
    # bc_umi_name_dict:  bc:[[umi1, name1, read1], [umi2, name2, read2]


print('Updating hammed adt ...')

bc_umi_name_dict = {}
c = 0

start = time.time()

for bcumi, adtdict in concat_dict.items():
    umi_adt_count = []
    umi_name_count = []
    
    bc = bcumi[:16]
    umi = bcumi[16:]
    
    for adt in adtdict:
    
        adt_updated = adt_all[adt] # adt_all[adt] is after ham
        if adt_updated != 'no_match' and adt_updated != 'two_matches':
            # Store umi, samplename, and read counts
            umi_adt_count = [umi, names[adt_updated], adtdict[adt]]
            if bc not in bc_umi_name_dict:
                bc_umi_name_dict[bc] = [umi_adt_count]
            else:
                temp_list = bc_umi_name_dict[bc]
                temp_list.append(umi_adt_count)
                bc_umi_name_dict[bc] = temp_list
            

    c+=1
    if c%500000 == 0:
        print('Finished %d reads' %c)
        
print('--- %d seconds ---' % ((time.time() - start)))






# 2. Find true bc by hamming (compared to the barcode.tsv ref)

print('bc hamming ...')

# Make bc_forhamming_dict from bc_umi_name_dict (keys)
bc_forhamming_dict = {}

for key in list(bc_umi_name_dict.keys()):
    bc_forhamming_dict[key] = 0
    

print (np.version.version)
print (sys.version)

start = time.time()

bc10x = make_bc10x(input_dir, bc_filename)
bc10x = np.array(bc10x)

hamAssociated = bc_forhamming_dict 

n=0
for bc in bc_forhamming_dict:
    hamAssociated[bc] = handleOne(bc, bc10x)
    n+=1
    if n%10000==0:
        print ('Finished %s reads' %n)

bc_meta_dict = {}
for h in hamAssociated:
    bc_meta_dict[h] = hamAssociated[h] # bc read = bc corrected
        
print('--- %d minutes ---' % ((time.time() - start)/60))


bc_corrected_dict = {}

for bc in bc_meta_dict:
    bc_updated = bc_meta_dict[bc] # bc after hamming
    
    if bc_updated != 'Remove_Minham' and bc_updated != 'Remove_multham1':
        if bc_updated not in bc_corrected_dict:
            bc_corrected_dict[bc_updated] = bc_umi_name_dict[bc]
        else:
            temp = bc_corrected_dict[bc_updated]
            temp.extend(bc_umi_name_dict[bc])
            bc_corrected_dict[bc_updated] = temp



# 3. UMI deconvolution with UMI clusterer #
# Reusing code for ecb here. RG = samplename

print('\numi deconvolution ...')

single_cleaned_dict = {} # to store {bc:[[umi, rg, count], [umi, rg, count], ...], ...}
    # only single rg is 'cleaned' in this dict
    
bc_name_assoc_dict = {} # to store confirmed bc:name pair

forprinting = {}

c = 0

if 'AACGGGAGTTCACGAT' in bc_corrected_dict:
    print(bc_corrected_dict['AACGGGAGTTCACGAT'])

for bc in bc_corrected_dict:
    umi_count_dict = {} # to count all umi within the same cell/droplet
    umi_rg_dict = {} # to store rg associated with the umi 
    umi_rg_count = {} # to store rg count associated with the umi
    rg_count_dict = {} # to keep track of different rg in the same cell (same bc)
    rg_umi_count = {}
    
    for entry in bc_corrected_dict[bc]:
        umi = entry[0]
        rg = entry[1]
        rg_count = entry[2]
        if umi not in umi_count_dict:
            umi_count_dict[umi] = 1
            umi_rg_dict[umi] = rg
            umi_rg_count[umi] = int(rg_count)
        else:
            if rg == umi_rg_dict[umi]:
                umi_rg_count[umi] += int(rg_count)

    byte_dict = {}  
    for seq in umi_count_dict:
        byte_dict[seq.encode()] = umi_count_dict[seq]
    
    seq_clustered = clusterer(byte_dict, 1) # return a list of clusters (3D list)
    
    # Go through every single cluster in the list 
    # We want to keep ALL clusters
    # First entry in each cluster is the 'true' umi
    for cluster in seq_clustered:
        rg_final_dict = {} # record rg:rg_count for each umi cluster
        umi_count = 0
        # Decode from byte to str
        cluster_string = [item.decode() for item in cluster]
        # true umi == first entry in each cluster 
        true_umi = cluster_string[0]
        for entry in cluster_string:
            if umi_rg_dict[entry] not in rg_final_dict: # rg associated with the umi. rg stored as int
                rg_final_dict[umi_rg_dict[entry]] = umi_rg_count[entry] # rg:rg_count
                rg_umi_count[umi_rg_dict[entry]] = umi_count_dict[entry] # rg:umi_count
            else:
                rg_final_dict[umi_rg_dict[entry]] += umi_rg_count[entry]
        
        # Sort rg_final_dict
        rg_sorted = sorted(rg_final_dict.items(), key=operator.itemgetter(1), reverse=True)
        
        # Keep the rg with the highest count: one umi should only have 1 rg
        if len(rg_sorted) == 1 or rg_sorted[0][1] > rg_sorted[1][1]:
            
            # track rg and count in the same cell
            if rg_sorted[0][0] not in rg_count_dict: 
                rg_count_dict[rg_sorted[0][0]] = rg_umi_count[rg_sorted[0][0]] # rg:umi_count
            else:
                rg_count_dict[rg_sorted[0][0]] += rg_umi_count[rg_sorted[0][0]]
               
            # store bc:[[umi, rg, rg_count], [umi2, rg, rg_count], ...]
            if bc not in single_cleaned_dict:
                single_cleaned_dict[bc] = [[true_umi, rg_sorted[0][0], rg_sorted[0][1]]]
            else:
                temp = single_cleaned_dict[bc]
                temp.append([true_umi, rg_sorted[0][0], rg_sorted[0][1]])
                single_cleaned_dict[bc] = temp  


# 4. Name cleanup: using bc_rg_set to track different rg for the same bc
    # In theory, there should only be 1 name (adt) associated with 1 bc

    ### Use umi count instead of reads
    
    total_count = 0 # umi count within a BC. Later, discard occurance < 2     
    true_rg_dict = {} # keep track of true rg count
    
    if len(rg_count_dict) == 1:
        for rg in rg_count_dict:
            total_count += rg_count_dict[rg] 
            true_rg_dict[rg] = rg_count_dict[rg]
        
    elif len(rg_count_dict) > 1: 
        total_count = sum(list(rg_count_dict.values()))

        # Set threshold = rg_count required to be considered a 'true rg'
        # true_name_threshold = 0.3 (Change at ## Variables ##)
        threshold = total_count/len(rg_count_dict)*true_name_threshold 
        #print(threshold)
        if total_count > rg_cleanup_thres: # rg_cleanup_thres = 2 (Change at ## Variables ##)
            for rg in rg_count_dict:
                if rg_count_dict[rg] > threshold and rg_count_dict[rg] >= umilim: # lousy reads will be removed
                    true_rg_dict[rg] = rg_count_dict[rg]

    
    # Filter: discard entry with less than total 20 umi count 
    if total_count < name_totalcount_lim:  # name_totalcount_lim = 20 
        if bc in single_cleaned_dict:
            single_cleaned_dict.pop(bc)
        bc_name_assoc_dict[bc] = 'NA'
        # store for printing: BC10X, Single_status, Total_Counts, Read_group, RG_Counts, Sample_Origin
        forprinting[bc] = [len(true_rg_dict), total_count, 'NA', 'NA', 'NA']
        
    else: # passed threshold
            
        if len(true_rg_dict) == 1:
            
            true_rg_local = list(true_rg_dict.keys())[0] # true rg as a string
            
            # Update the rg in single_cleaned_dict and bc_rg_assoc_dict 
            entry_list = single_cleaned_dict[bc] #[[true_umi, rg, rg_count], ...]
            entry_true_rg = [] # to store entries with true rg
            for entry2 in entry_list: # [true_umi, rg, rg_count]
                if entry2[1] == true_rg_local: # if rg_sorted[0][0] is a true rg
                    entry_true_rg.append(entry2) # keep the entry
            single_cleaned_dict[bc] = entry_true_rg # now only 1 rg is associated to each bc
            bc_name_assoc_dict[bc] = true_rg_local
            
            # store for printing: BC10X, Single_status, Total_Counts, Read_group, RG_Counts, Sample_Origin
            name_info = list(true_rg_dict.items())[0]
            forprinting[bc] = [len(true_rg_dict), total_count, name_info[0], name_info[1], name_info[0]]
        

        
        # more than one true name
            # Check if any name dominate (exceed 70% of total count) -> assign as sample_origin
            # else: multiplet
        elif len(true_rg_dict) > 1:
            name_linked = []
            count_linked = []
            sample_origin = 'Multiplet'
            
            # sort true_rg_dict in descending order by value (count)
            name_sorted = sorted(true_rg_dict.items(), key=operator.itemgetter(1), reverse=True)
            # major_name_percent = 0.7 (Change at ## Variables ##)
            if name_sorted[0][1] > total_count*major_name_percent: 
                sample_origin = name_sorted[0][0]
            
            for nametuple in name_sorted:
                name_linked.append(nametuple[0])
                count_linked.append(str(nametuple[1]))
            
            # store for printing: BC10X, Single_status, Total_Counts, Read_group, RG_Counts, Sample_Origin
            forprinting[bc] = [len(true_rg_dict), total_count, ';'.join(name_linked),
                               ';'.join([str(x) for x in count_linked]), sample_origin]
            
            # store sample_origin in bc_name_assoc_dict
            bc_name_assoc_dict[bc] = sample_origin
            
            if sample_origin == 'Multiplet':
                single_cleaned_dict.pop(bc)
    c+=1
    if c%5000 == 0:
        print('Finished %d reads' %c)



# Sort forprinting dict by key (bc)
forprinting_sorted = {key:forprinting[key] for key in sorted(forprinting.keys())}



## Intermediate Output ##
# bc_adt_count_outfile: Hash_correlation.txt
# equvalent to Hash_correlation_condensation_clean_min20all_more70p.txt

header = ['BC10X','Singlet_status','Total_counts','ADT_group','ADT_group_counts','Sample_Origin']

with open(os.path.join(interpath, bc_adt_count_outfile), 'w') as outf:
    print('\t'.join(header), file=outf)
    for bc, entry in forprinting_sorted.items():
        print(bc+'\t'+'\t'.join([str(x) for x in entry]), file=outf)
        

#### Output ####
# seurat_outfile: 'Hash_forSeurat.csv'
# equivalent to Hash_barcodes_for_Seurat_min20all_more70p
    # Use bc_name_assoc_dict. bc: name
    # Include all bc from the whitelist (bcwhitelist)

print('creating output %s' %seurat_outfile)

bcwhitelist = []
seuratdict = {}
name_list = []

bc_file = open(os.path.join(input_dir, bc_filename), 'rt')
for line in bc_file:
    bcwhitelist.append(line.strip('\n').rstrip('-1'))  
    
for bc, name in bc_name_assoc_dict.items():
    if name != 'Multiplet' and name != 'NA':
        RGname = 'RG_' + name
        seuratdict[bc] = RGname # bc in meta_dict_update should be unique
        if RGname not in name_list:
            name_list.append(RGname)
    elif name == 'Multiplet':
        seuratdict[bc] = 'RG_Multiplet'
    elif name == 'NA':
        seuratdict[bc] = 'RG_NA'    
        
        
name_list.append('RG_Multiplet')
name_list.append('RG_NA')


with open(os.path.join(output_dir, seurat_outfile), 'w') as outf:
    w = csv.writer(outf)
    w.writerow(bcwhitelist)
    for name2 in name_list:
        temp_row = []
        temp_row.append(name2)
        for bc2 in bcwhitelist:
            if bc2 in seuratdict:
                if seuratdict[bc2] == name2:
                    temp_row.append(1000)
                else:
                    temp_row.append(0)
            else: 
                temp_row.append(0)

        w.writerow(temp_row)



print('--- %d minutes ---'%((time.time()-BeginningOfTime)/60))


