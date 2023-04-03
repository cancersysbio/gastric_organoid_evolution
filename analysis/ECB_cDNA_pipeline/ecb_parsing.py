#!/usr/bin/python -u
# coding: utf-8


## GENERAL ABBREVIATIONS ##
# ecb = ExtraCellular Barcode, which is the DNA barcode
# bc = 10X BarCode
# umi = unique molecular identifier




import os
import sys
import gzip
from umi_tools.network import UMIClusterer
clusterer = UMIClusterer(cluster_method="directional")
from ecb_functions import * # uncomment for bash
import time
import copy
import operator
import numpy as np
import csv
from collections import Counter
from datetime import date


##----------------------Variables------------------------------##
input_dir = sys.argv[1]
fastq_name = sys.argv[2]
ecb_master = sys.argv[3]

## Set Paths ##
# Below paths are auto #
output_dir = input_dir + 'ecb_outs/'
if not os.path.exists(output_dir):
    os.mkdir(output_dir) 
interpath = output_dir + 'intermediate_output/'
if not os.path.exists(interpath):
    os.mkdir(interpath)
   
## fastq parsing parameters ##
ecb_handle = 'CTGCAGTCTGAGTCTGACAG'
hamming_distance = 2
majority_percent = 50 # Used 50 in gastric tumor evolution

## thresholds and ratio for filtering ##
umi_min_read = 2 # used in step 3, default 2, (used 5 for seq 22 and 23, high read depth)
rg_cleanup_thres = 2 # used in step 4
umilim = 2 # used in step 4, default 2 (used 5 for seq 22 and 23, high read depth)
truerg_thres_low = 0.3 # used in step 4
truerg_thres_high = 0.1 # used in step 4 (used 0.1 in gastric tumor evolution)
max_rg = 4 # used in step 4
rg_multicount_thres = 10 # step 5
multi_single_rg_ratio = 0.15 # step 5, default 0.15 (0.5 for D3C2, 1 for D2C5R1T8_LATE)

## ecb and bc reference list ##
#ecb_ref = 'ECB_groups_master_list.txt'
ecb_ref = ecb_master + '.txt'
bc_filename = 'barcodes.tsv'

#### Intermediate Output ####  
ecb_outfile = 'ECB_ham2_clean_RG.txt'
lost_outfile = 'umi_lost_after_deconvolution.txt'
forSeurat_outfile = 'forSeurat.csv'
BC_pairs_outfile = 'BC_pairs.txt'
true_pairs_outfile = 'BC_pairs_true.txt'

#### Outputs ####
step4meta_outfile = 'BC_RG_AssociatedUMICount.txt'
newrg_outfile = 'BC_RG_AssociatedUMICount_updated.txt'
concat_rg_outfile = 'BC_pairs_true_multi.txt'
forSeurat_updated_outfile = 'forSeurat_updated.csv'
ecb_ref_updated_outfile = 'ECB_groups_master_list_updated.txt'
    
##------------------------------------------------------------##


BeginningOfTime = time.time() # Track time for the whole run


## Parse ecb fastq ##
concat_dict = {}

if ',' not in fastq_name: # single run
    print('Parsing ecb fastq...')
    concat_dict = extract_fastQ(input_dir, fastq_name, ecb_handle, hamming_distance) #{bc_umi:{ecb1:count, ecb2:count ...}}
    print('Parsing ecb done.')
else: # multiple runs
    c = 0
    fastq_name = fastq_name.split(',')
    for fastq in fastq_name:
        print('Parsing ecb %s...' %fastq)
        c += 1
        temp_dict = {}
        temp_dict = extract_fastQ(input_dir, fastq, ecb_handle, hamming_distance)
        print(len(temp_dict))
        if c == 1: # first pair of R1R2
            concat_dict = {**temp_dict}
           # concat_dict = copy.deepcopy(temp_dict) 
        elif c > 1: # second pair and so on
            for key in temp_dict:
                if key not in concat_dict:
                    concat_dict[key] = temp_dict[key]
                else:
                    subdict = concat_dict[key] # {ecb:read}
                    subdict = Counter(subdict) + Counter(temp_dict[key])
                    concat_dict[key] = subdict
        print('%s done' %fastq)
len(concat_dict)




## Deconvolution ##
# 1. Find major ecb by UMIclusterer
    # Combine count of all valid ecb and associate to the major ecb
    # Assign RG index

print('ecb deconvolution...')

# Parse ecb_ref to create a ecb:readgroup index dict
ecb_ref_dict = {}

with open(os.path.join(input_dir, ecb_ref), 'rt') as ecb_f:
    for line in ecb_f:
        ecb, rg = line.strip().split('\t')
        ecb_ref_dict[ecb] = rg

start_time = time.time()

total = 0
rgmatch = 0
norg = 0
lowsum = 0
lowmajor = 0

bc_umi_rg_dict = {}
ecb_print_list = [] # storing info for printing. List not dict because bc repeats
ecb_lost_list = []

for bc_umi in concat_dict:
    ecb_rg_dict = {}
    total += 1
    
    bc = bc_umi[:16] 
    umi = bc_umi[16:]
    
    ecb_dict = concat_dict[bc_umi]
    
    # count all ecb
    ecb_sum = 0
    for k in ecb_dict:
        ecb_sum += ecb_dict[k]

    # Sort ecb_dict
    ecb_sorted_dict = {ecb: count for ecb,count in sorted(ecb_dict.items(), key=lambda item: item[1], reverse=True)}

    # Assign rg to every single ecb (absolute matching)
    for ecb in ecb_sorted_dict:
        if ecb in ecb_ref_dict:
            ecb_rg_dict[ecb] = ecb_ref_dict[ecb] # ecb:rg    
        else:
            ecb_rg_dict[ecb] = -1
        
    major_rg = -1
    # Assign major ecb and major rg
    for ecb2 in ecb_sorted_dict:
        if major_rg == -1: # rg not assigned yet
            if ecb_rg_dict[ecb2] != -1: 
                ## These should only happen once per bc ##
                rgmatch += 1
                major_ecb = ecb2
                major_rg = ecb_rg_dict[ecb2]
                major_count = ecb_sorted_dict[ecb2]
    
        else: # if major_rg has already been assigned
            if ecb_rg_dict[ecb2] == major_rg: # if same rg
                major_count += ecb_sorted_dict[ecb2] # combine count
    
    
    if major_rg != -1: # if a rg has been assigned
        if major_count >= (ecb_sum*majority_percent/100): # and its count is greater than or equal to 50% of the total ecb count
            # Keep and store info for printing
            if ecb_sum > 1:
                # store info for printing intermediate output
                ecb_print_list.append([bc, umi, major_ecb, major_count, ecb_sum, major_rg])

                # store info for the next step
                if bc not in bc_umi_rg_dict:
                    bc_umi_rg_dict[bc] = [[umi, major_rg, major_count]]
                else:
                    temp_list = bc_umi_rg_dict[bc]
                    temp_list.append([umi, major_rg, major_count])
                    bc_umi_rg_dict[bc] = temp_list

            else: # if ecb_sum is 1
                lowsum += 1
                for ecb_failed in ecb_sorted_dict:
                    ecb_lost_list.append([bc, umi, ecb_failed, 
                                          ecb_sorted_dict[ecb_failed], 'ecb_sum_1'])
        else: # if major_count is too low
            lowmajor += 1

    elif major_rg == -1: # if no rg assigned
        norg += 1
        for ecb_failed in ecb_sorted_dict:
            ecb_lost_list.append([bc, umi, ecb_failed, 
                                  ecb_sorted_dict[ecb_failed], 'no_rg'])
     
            
    if total%100000 ==0:
         print('Finished %d reads' % total) 
            
print('--- %d seconds ---' % (time.time() - start_time))

print('Started with %s entries. %s matched with a rg.\n'
      '%s dropped due to rg count lower than 50 percent of total count.\n'
      '%s dropped due to no rg match.\n'
      '%s dropped due to major ecb less than 2.\n'%(total, rgmatch, lowmajor, norg, lowsum))


#### Intermediate Output ####

# Output: Equivalent to Valid_cellBC_ECB_ham2_clean_RG.txt
 # Use info stored in ebc_print_list
 # 6 Columns: bc, umi, ecb, major_ecb_count, all_ecb_count, RG


# ecb_outfile = 'ECB_ham2_clean_RG.txt'
with open(os.path.join(interpath, ecb_outfile), 'w') as outf:
    header = ['bc', 'umi', 'ecb', 'major_ecb_count', 'all_ecb_count', 'RG']
    print('\t'.join(header), file=outf)
    for sublist in ecb_print_list:
        print('\t'.join([str(x) for x in sublist]), file=outf)
      
    
# lost_outfile = 'umi_lost_after_deconvolution.txt'
with open(os.path.join(interpath, lost_outfile), 'w') as lostf:
    header2 = ['bc', 'umi', 'ecb', 'ecb_count', 'fail_reason']
    print('\t'.join(header2), file=lostf)
    for sublt in ecb_lost_list:
        print('\t'.join([str(y) for y in sublt]), file=lostf)








# 2. Find true bc by hamming (compared to the barcode.tsv ref) 

print('10xbc hamming...')

bc_forhamming_dict = {}
for bc in bc_umi_rg_dict:
    bc_forhamming_dict[bc] = 0


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
        
print('--- %d minutes ---' % ((time.time() - start_time)/60))


bc_corrected_dict = {}

for bc in bc_meta_dict:
    if bc_meta_dict[bc] != 'Remove_Minham' and bc_meta_dict[bc] != 'Remove_multham1':
        if bc_meta_dict[bc] not in bc_corrected_dict:
            bc_corrected_dict[bc_meta_dict[bc]] = bc_umi_rg_dict[bc]
        else:
            temp = bc_corrected_dict[bc_meta_dict[bc]]
            temp.extend(bc_umi_rg_dict[bc])
            bc_corrected_dict[bc_meta_dict[bc]] = temp



# 3. UMI deconvolution with UMI clusterer #
print('\numi deconvolution...\n')

single_cleaned_dict = {} # to store {bc:[[umi, rg, count], [umi, rg, count], ...], ...}
    # only single rg is 'cleaned' in this dict
    
bc_rg_assoc_dict = {} # to store confirmed bc:rg pair

### these 3 dicts store info for step 5: true rg pair ###
single_rg_dict = {} # cell count of all significant single rg
multi_rg_dict = {} # to store all potential multi rg and their cell count
double_rg_dict = {} # to store all splitted up multi rg
bc_multi_rg_dict = {} # to store bc associated with the potential multi rg

forprinting = {}

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
            # rg associated with the umi. rg stored as int
            if int(umi_rg_dict[entry]) not in rg_final_dict:
                rg_final_dict[int(umi_rg_dict[entry])] = umi_rg_count[entry] # rg:rg_count
                rg_umi_count[int(umi_rg_dict[entry])] = umi_count_dict[entry] # rg:umi_count
            else:
                rg_final_dict[int(umi_rg_dict[entry])] += umi_rg_count[entry]
        
        # Sort rg_final_dict
        rg_sorted = sorted(rg_final_dict.items(), key=operator.itemgetter(1), reverse=True)
        

        # Keep the rg with the highest count: one umi should only have 1 rg
        if rg_sorted[0][1] >= umi_min_read: # keep only if there are at least 2 reads support the umi
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


# 4. RG cleanup: using bc_rg_set to track different rg for the same bc
    # In theory, there should only be 1 rg associated with 1 bc
        # which is one cell (rg) in on droplet (bc)
    # If 1 bc has more than 1 rg, usually it infers two cells trapped in the same droplet
        # or 1 cell was infected with more than 1 ecb

    ### Use umi count instead of reads
    
    total_count = 0 # umi count within a BC. Later, discard occurance < 2     
    true_rg_dict = {} # keep track of true rg count
    
    if len(rg_count_dict) == 1:
        for rg in rg_count_dict:
            total_count += rg_count_dict[rg] 
            true_rg_dict[rg] = rg_count_dict[rg]
        
    elif len(rg_count_dict) > 1: 
        for rg in rg_count_dict:
            total_count += rg_count_dict[rg] 

        # Set threshold = rg_count required to be considered a 'true rg'
            # For example: {'3': 53, '0': 68, '2': 3, '1': 3}
            # '3' and '0' are true rg but not '2' and '1'
        # truerg_thres_low = 0.3. (change at ### Variables ###)
        threshold = total_count/len(rg_count_dict)*truerg_thres_low
        # truerg_thres_high = 0.2 (change at ### Variables ###)
        threshold_highcount = total_count/len(rg_count_dict)*truerg_thres_high 
        if total_count > rg_cleanup_thres: # rg_cleanup_thres = 2 (change at ### Variables ###)
            for rg in rg_count_dict:
                if rg_count_dict[rg] > threshold and rg_count_dict[rg] >= umilim: # lousy reads will be removed
                    true_rg_dict[rg] = rg_count_dict[rg]
                elif (total_count > 100 and rg_count_dict[rg] >= threshold_highcount 
                      and rg_count_dict[rg] >= umilim):
                    true_rg_dict[rg] = rg_count_dict[rg]


    # Filter: discard entry with less than total 2 umi count (change at cell 2)
    if total_count < rg_cleanup_thres: 
        if bc in single_cleaned_dict:
            single_cleaned_dict.pop(bc)
#    elif len(true_rg_dict) > 4:  # Or remove lousy reads
#        if bc in single_cleaned_dict:
#            single_cleaned_dict.pop(bc)
    elif not true_rg_dict: # if true_rg_dict is empty
        if bc in single_cleaned_dict:
            single_cleaned_dict.pop(bc)
        
    else: # passed threshold
        
        # store for printing
        forprinting[bc] = true_rg_dict
        
        if len(true_rg_dict) == 1:
            true_rg_local = list(true_rg_dict.keys())[0] # true rg as a string
        
            if true_rg_local in single_rg_dict: # count single rg for step 5 (true pair testing)
                single_rg_dict[true_rg_local] += 1 # cell count
            else:
                single_rg_dict[true_rg_local] = 1 # cell count

            # Update the rg in single_cleaned_dict and bc_rg_assoc_dict 
            entry_list = single_cleaned_dict[bc] #[[true_umi, rg, rg_count], ...]
            entry_true_rg = [] # to store entries with true rg
            for entry2 in entry_list: # [true_umi, rg, rg_count]
                if entry2[1] == true_rg_local: # if rg_sorted[0][0] is a true rg
                    entry_true_rg.append(entry2) # keep the entry
            single_cleaned_dict[bc] = entry_true_rg # now only 1 rg is associated to each bc
            bc_rg_assoc_dict[bc] = true_rg_local
        
        
        # more than one true rg
            # store rg1_rg2_rg3: total count
        elif len(true_rg_dict) > 1 and len(true_rg_dict) <= max_rg:
            # use list so the rgs can be sorted. Tuple for storing as key in dict
            rg_multi = tuple(sorted(list(true_rg_dict.keys()))) 

            if rg_multi in multi_rg_dict:
                multi_rg_dict[rg_multi] += 1 # cell count
            else:
                multi_rg_dict[rg_multi] = 1 # cell count
            
            tempdouble = []
            for rg1 in rg_multi:
                for rg2 in rg_multi:
                    if rg1 != rg2:
                        temppair = tuple(sorted([rg1, rg2]))
                        if temppair not in tempdouble:
                            tempdouble.append(temppair)
                            
            for pair in tempdouble:
                if pair not in double_rg_dict:
                    double_rg_dict[pair] = 1
                else:
                    double_rg_dict[pair] += 1

            # Store this bc + multi combo as a reference later
            bc_multi_rg_dict[bc] = rg_multi



#### Output ####
# step4meta_outfile: BC_RG_AssociatedUMICount.txt
# equivalent to Valid_cellBC_ECB_ham2_clean_RG_test2_collapse3_umilim_1.txt
# bc, rg, umi_count, rg_count
    # make an output with: bc (single_cleaned_dict), rg (single_cleaned_dict+bc_multi_rg_dict), 
    #+ rg-umi.count(single_cleaned_dict) , rg_count (single_cleaned_dict)
    #+ create a umi rg dict with each bc
    #+ use forprinting as a ref for multi
    
print('creating output', step4meta_outfile)

forprint1 = []
meta_dict = {} # store all printed info
header = ['bc', 'rg', 'umi_count', 'read_count']

for bc in single_cleaned_dict:
    rg_umi_print = {}
    rg_reads = {}
    long_entry = single_cleaned_dict[bc]
    if forprinting[bc]:
    
        if bc not in bc_multi_rg_dict:

            for entry in long_entry:
                if entry[1] not in rg_umi_print: # rg
                    rg_umi_print[entry[1]] = 1 # umi count
                else:
                    rg_umi_print[entry[1]] += 1
                if entry[1] not in rg_reads:
                    rg_reads[entry[1]] = entry[2] # read count
                else:
                    rg_reads[entry[1]] += entry[2] 

        else:
            for m_rg in bc_multi_rg_dict[bc]:
                for entry in long_entry:
                    if m_rg == entry[1]:
                        if m_rg not in rg_umi_print: # rg
                            rg_umi_print[m_rg] = 1 # umi count
                        else:
                            rg_umi_print[m_rg] += 1
                        if m_rg not in rg_reads:
                            rg_reads[m_rg] = entry[2] # read count
                        else:
                            rg_reads[m_rg] += entry[2]


        if len(rg_umi_print) == 1:
            rgprint = list(rg_umi_print.keys())[0]
            forprint1.append('\t'.join([bc, str(rgprint), str(rg_umi_print[rgprint]), str(rg_reads[rgprint])]))

        else:
            forprint1.append('\t'.join([bc, '_'.join([str(x) for x in list(rg_umi_print.keys())]), 
                    '_'.join([str(y) for y in list(rg_umi_print.values())]), 
                    '_'.join([str(z) for z in list(rg_reads.values())])]))
            
with open(os.path.join(output_dir, step4meta_outfile), 'w') as outf:
    print('\t'.join(header), file=outf)
    for line in forprint1:
        print(line, file=outf)
        meta_dict[line.strip().split('\t')[0]] = line.split('\t')[1:]
        


#### Intermediate Output ####

# forSeurat_outfile: equivalent to Valid_cellBC_ECB_ham2_clean_RG_test2_collapse3_umilim_1_forSeurat.csv
    # Use bc_rg_assoc_dict and bc_multi_rg_dict
    # any multi RG will be 'duplet' now
    # Include all bc from the whitelist (bc_filename)
    
print('creating output', forSeurat_outfile)
    
bcwhitelist = []
seuratdict = {}
rg_list = []

bc_file = open(os.path.join(input_dir, bc_filename), 'rt')
for line in bc_file:
    bcwhitelist.append(line.strip('\n').rstrip('-1'))

for bc in bc_rg_assoc_dict:
    if bc not in bc_multi_rg_dict:
        rg_n = '_'.join(['RG', str(bc_rg_assoc_dict[bc])])
        seuratdict[bc] = rg_n
        if rg_n not in rg_list:
            rg_list.append(rg_n)
    else:
        seuratdict[bc] = 'RG_duplet'

rg_list.append('RG_duplet')

with open(os.path.join(interpath, forSeurat_outfile), 'w') as outf:
    w = csv.writer(outf)
    w.writerow(bcwhitelist)
    for rg_n2 in rg_list:
        temp_row = []
        temp_row.append(rg_n2)
        for bc2 in bcwhitelist:
            if bc2 in seuratdict:
                if seuratdict[bc2] == rg_n2:
                    temp_row.append(1000)
                else:
                    temp_row.append(0)
            else: 
                temp_row.append(0)
        w.writerow(temp_row)




#### Intermediate Output ####

# BC_pairs_outfile - BC_pairs.txt
# equivalent to Valid_cellBC_ECB_ham2_clean_RG_test2_collapse3_umilim_1_multi_BC_pairs.txt
    # cell count of and single and double rg
    #+ use single_rg_dict and double_rg_dict

header = ['rg', 'number_of_cells_associated']
    
with open(os.path.join(interpath, BC_pairs_outfile), 'w') as outf:
    print('\t'.join(header), file=outf)
    for s_rg in single_rg_dict:
        print(s_rg, '\t', single_rg_dict[s_rg], file=outf)
    for multirg in double_rg_dict:
        multi_rg = '_'.join([str(x) for x in list(multirg)])
        print(multi_rg, '\t', double_rg_dict[multirg], file=outf)



# 5. Check if the multi_rg are true pair
    # test if the ratio of multi_rg cell count to single rg cell count is high enough 
        # change multi_single_rg_ratio at cell 2. Default multi_single_rg_ratio = 0.15 
    # combine any double with overlapping rg
        
multi_rg_meta = {} # store info and outcome as a reference
multi_cleaned_dict = {} # store completely clean {bc:[[umi, rg, count], [umi, rg, count], ...], ...}
true_multi_rg_dict = {}

# copy single rg data to multi_cleaned_dict
for bc in single_cleaned_dict:
    if bc not in bc_multi_rg_dict:
        multi_cleaned_dict[bc] = single_cleaned_dict[bc]
    
for multi_rg in double_rg_dict: # multi_rg are tuples
    multi_sum = double_rg_dict[multi_rg] # number of cells supporting this multi rg
    rg1 = multi_rg[0]
    rg2 = multi_rg[1]
    rg1_count = 0
    rg2_count = 0
    if rg1 in single_rg_dict:
        rg1_count = single_rg_dict[rg1]    
    if rg2 in single_rg_dict:
        rg2_count = single_rg_dict[rg2]
    
    single_sum = rg1_count + rg2_count
    go = 'no' # a switch for true multi rg
    
    # test if multisum/singlesum >= and multisum >= 10
    if single_sum > 0: # it is possible that rg in multi_rg do not exist in single
        # Default multi_single_rg_ratio = 0.15, rg_multicount_thres = 10. Change at cell 2
        if multi_sum/single_sum >= multi_single_rg_ratio and multi_sum > rg_multicount_thres:
            multi_rg_meta[multi_rg] = [str(multi_sum), str(single_sum), str(rg1_count),
                                       str(rg2_count), str(multi_sum/single_sum), 'TRUE'] 
                                        # store as a reference
            # store valid multi rg and their cell count
            if multi_rg not in true_multi_rg_dict:
                true_multi_rg_dict[multi_rg] = multi_sum
            else:
                true_multi_rg_dict[multi_rg] += multi_sum
            
    elif single_sum == 0 and multi_sum > rg_multicount_thres:
        multi_rg_meta[multi_rg] = [str(multi_sum), 'no single rg']
        if multi_rg not in true_multi_rg_dict:
            true_multi_rg_dict[multi_rg] = multi_sum
        else:
            true_multi_rg_dict[multi_rg] += multi_sum
            



# 6. Combine overlapping true multi_RG
    # First, create an union_rg list based on true_multi_rg_dict
    
# sort true_multi_rg_dict
true_multi_rg_dict = {key:true_multi_rg_dict[key] for key in sorted(true_multi_rg_dict)}

# store all multi_rg in a list
mrg_list = list(true_multi_rg_dict.keys())

# convert each multi_rg from tuple to set
mrg_sets = [set(mrg) for mrg in mrg_list if mrg]
merge = True
while merge:
    merge = False
    combined_list = []
    while mrg_sets:
        overlap, rest = mrg_sets[0], mrg_sets[1:]
        mrg_sets = []
        for mrg2 in rest:
            if mrg2.isdisjoint(overlap):
                mrg_sets.append(mrg2)
            else:
                merge = True
                overlap = {*overlap, *mrg2}
        combined_list.append(overlap)
    mrg_sets = combined_list
union_rg_list = mrg_sets # stored




### Intermediate Output ###
    # true_pairs_outfile
    # BC_pairs_true.txt: equivalent to Valid_cellBC_ECB_ham2_clean_RG_test2_collapse3_umilim_1_multi_BC_pairs_true.tx
    # Use multi_rg_meta: rg1_rg2, doublesum, singlesum, rg1_singlecount, rg2_singlecount, ratio, 'TRUE'
    
header = ['rg', 'doublesum', 'singlesum', '1st_rg_sum', '2nd_rg_sum', 'ratio', 'truepair']
    
with open(os.path.join(interpath, true_pairs_outfile), 'w') as outf:
    print('\t'.join(header), file=outf)
    for mrg in multi_rg_meta:
        print('_'.join([str(x) for x in mrg]), '\t', '\t'.join([str(y) for y in multi_rg_meta[mrg]]), file=outf)



### Output ###
    # concat_rg_outfile
    # BC_pairs_true_multi.txt: equivalent to Valid_cellBC_ECB_ham2_clean_RG_test2_collapse3_umilim_1_multi_BC_pairs_true_multi 
    # union_rg_list
    # multi_rg \t lowest rg
    
print('creating output', concat_rg_outfile)

header = ['multi_rg', 'lowest_rg']

with open(os.path.join(output_dir, concat_rg_outfile), 'w') as outf:
    print('\t'.join(header), file=outf)
    for mrg in union_rg_list:
        mrg_sorted = sorted(mrg)
        print('_'.join([str(x) for x in mrg_sorted]), '\t', mrg_sorted[0], file=outf)
        


# 7. Collapse
# e.g. (14, 20): 57
# rg20 and rg14_rg20 will be renamed into rg14 but NOT rg14_rg20_rg22
# use true_multi_rg_dict_sorted, multi_rg_dict, bc_rg_assoc_dict -> multi_rg_collapsed_dict
# Then bc_multi_rg_dict?

rg_collapsed_dict = {}

for true_m in union_rg_list: # (14, 20)
    true_m = sorted(list(true_m)) # because sets are not iterable
#     print(true_m, min(true_m), true_m[0])
    # Go through all multi_rg
    for mrg in multi_rg_dict: # e.g. (0, 14, 20)
        match_count = 0
        for ori_rg in mrg: # 0, 14, 20
            if ori_rg in true_m:
                match_count += 1
        if match_count == len(mrg):
            rg_collapsed_dict[mrg] = min(true_m) # assign lowest rg
        else:
            if mrg not in rg_collapsed_dict:
                rg_collapsed_dict[mrg] = mrg # no change
        
    # Go through all single_rg
    for bc in bc_rg_assoc_dict:
        srg = bc_rg_assoc_dict[bc]
        if srg in true_m:
            rg_collapsed_dict[srg] = min(true_m) # assign lowest rg
        else: 
            if srg not in rg_collapsed_dict:
                rg_collapsed_dict[srg] = srg # no change



# copy ecb rg ref
ecb_paired_dict = {**ecb_ref_dict}



## Output ##
# ecb_ref_updated_outfile: ECB_groups_master_list_pairs.txt
# update ecb master list with updated rg after collaspe

for ecb, rg in ecb_ref_dict.items():
    for oldrg, newrg in rg_collapsed_dict.items():
        if type(oldrg) == int and int(rg) == oldrg: # single rg are stored as integers
            ecb_paired_dict[ecb] = str(newrg)

with open(os.path.join(output_dir, ecb_ref_updated_outfile), 'w') as outf:
    for key, item in ecb_paired_dict.items():
        print(key+'\t'+item, file=outf)




### Output ###
# newrg_outfile: BC_RG_AssociatedUMICount_update.txt
# equivalent to Valid_cellBC_ECB_ham2_clean_RG_test2_collapse3_umilim_1_pairs.txt
# bc, rg, umi_count, read_count, updated_rg
# use meta_dict and rg_collapsed_dict

print('creating output', newrg_outfile)

meta_dict_update = {}
header = ['bc', 'rg', 'umi_count', 'read_count', 'updated_rg']

for bc in meta_dict:
    rg, umicount, reads = meta_dict[bc]
    if '_' in rg:
        rg_int = tuple(int(x) for x in rg.split('_'))
    else:
        rg_int = int(rg)
        
    if rg_int in rg_collapsed_dict:
        new_rg = rg_collapsed_dict[rg_int]
        if type(new_rg) == tuple:
            new_rg_str = '_'.join([str(x) for x in new_rg])
        elif type(new_rg) == int:
            new_rg_str = str(new_rg)
        else:
            print('error')
        meta_dict_update[bc] = [rg, umicount, reads, new_rg_str]
    else:
        print(rg_int)

    if not rg_collapsed_dict: # if there are no true multi_rg at all
        meta_dict_update[bc] = [rg, umicount, reads, rg]

with open(os.path.join(output_dir, newrg_outfile), 'w') as outf:
    print('\t'.join(header), file=outf)
    for bc2 in meta_dict_update:
        print(bc2+'\t', '\t'.join(meta_dict_update[bc2]), file=outf)




#### Output ####

# forSeurat_updated_outfile: equivalent to Valid_cellBC_ECB_ham2_clean_RG_test2_collapse3_umilim_1_forSeurat.csv
    # Use meta_dict_update. bc:[rg, umi_count, read_count, new_rg]
    # Include all bc from the whitelist (bcwhitelist)

print('creating output', forSeurat_updated_outfile)
    
seurat_dict_updated = {}
rg_list2 = []
    
for bc in meta_dict_update:
    rg = meta_dict_update[bc][3]
    if '_' not in rg:
        rg_n = 'RG_' + rg
        seurat_dict_updated[bc] = rg_n # bc in meta_dict_update should be unique
        if rg_n not in rg_list2:
            rg_list2.append(rg_n)
    else:
        seurat_dict_updated[bc] = 'RG_duplet'
rg_list2.append('RG_duplet')

with open(os.path.join(output_dir, forSeurat_updated_outfile), 'w') as outf:
    w = csv.writer(outf)
    w.writerow(bcwhitelist)
    for rg_n2 in rg_list2:
        temp_row = []
        temp_row.append(rg_n2)
        for bc2 in bcwhitelist:
            if bc2 in seurat_dict_updated:
                if seurat_dict_updated[bc2] == rg_n2:
                    temp_row.append(1000)
                else:
                    temp_row.append(0)
            else: 
                temp_row.append(0)
        w.writerow(temp_row)





EndOfTime = time.time() - BeginningOfTime
if EndOfTime < 60:
    print('\nTotal time: %d seconds' %EndOfTime)
else:
    print('\nTotal time: %d minutes' %(EndOfTime/60))
