#!/usr/bin/env python
# coding: utf-8


import sys
import os



## Paths ##

ecb_dir = sys.argv[1]
hash_dir = sys.argv[2]
ref_dir = sys.argv[3]

out_dir = ecb_dir + 'combined_outs/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir) 
    
## Input files names ##
bcwhitelist = 'barcodes.tsv'
color_file = 'barcode_colors_hexcode_RG_match.txt'
bc_rg_file = 'BC_RG_AssociatedUMICount_updated.txt'
bc_adt_file = 'Hash_correlation.txt'

## Output files names ##
outfile = 'Cell_Metadata_updated.txt'
outfile_min10 = 'Cell_Metadata_min10cells_perECBrg_updated.txt'
outfile_full = 'Cell_Metadata_full_updated.txt'
outfile_summary = 'Cell_Metadata_summary_updated.txt'
outfile_ecbhash = 'Cell_Metadata_ecbs_per_hashSample_updated.txt'


# Create a BC ref list from the whitelist
bc_ref_list = []

with open(os.path.join(ref_dir, bcwhitelist), 'rt') as bfile:
    for line in bfile:
        bc = line.strip()[:-2]
        bc_ref_list.append(bc)


# Color reference
    # rg:color {dict}
rg_color = {}
    
with open(os.path.join(ref_dir, color_file), 'rt') as cfile:
    for line in cfile:
        rg, color = line.strip().split('\t')
        rg_color[rg] = color
  

# ECB
    # 1. RG cell count with 10 minimum cell count {dict}
    # 2. BC:RG dict {dict}
    # 3. BC list [list]
    # 4. RG:BC {dict}
# use BC_RG_AssociatedUMICount_updated.txt: bc, rg, umi_count, read_count, updated_rg

RG_cellcount = {}
RG_cellcount_min10 = {}
BC_RG = {}
bc_read_list = []
multiplet_count = 0
RG_BC = {}

with open(os.path.join(ecb_dir, bc_rg_file), 'rt') as efile:
    efile.readline()
    for line in efile:
        bc, rg, umi_count, read_count, updated_rg  = line.strip().split('\t')
        # RG cell count
        if updated_rg not in RG_cellcount:
            RG_cellcount[updated_rg] = 1
        else:
            RG_cellcount[updated_rg] += 1
            
        # BC:RG
        if '_' in updated_rg:
            multiplet_count += 1
            if 'Multiplet' not in RG_BC:
                RG_BC['Multiplet'] = [bc]
            else:
                templist = RG_BC['Multiplet']
                templist.append(bc)
                RG_BC['Multiplet'] = templist
        else:
            if bc not in BC_RG:
                BC_RG[bc] = updated_rg
            else:
                print('collision!')
            if updated_rg not in RG_BC:
                RG_BC[updated_rg] = [bc]
            else:
                templist = RG_BC[updated_rg]
                templist.append(bc)
                RG_BC[updated_rg] = templist
        
        # bc_read_list
        bc_read_list.append(bc)
            

for RG, count in RG_cellcount.items():
    if count > 10:
        RG_cellcount_min10[RG] = count
            
            
# ADT
    # bc:hash {dict}
bc_hash = {}

with open(os.path.join(hash_dir, bc_adt_file), 'rt') as afile:
    afile.readline() # skip header
    for line in afile:
        bc = line.rstrip().split('\t')[0]
        hash_tag = line.strip().split('\t')[5]
        bc_hash[bc] = hash_tag


# Combine #
    # outfile = 'Cell_Metadata.txt'
    # outfile_min10 = 'Cell_Metadata_min10cells_perECBrg.txt'
    # outfile_full = 'Cell_Metadata_full.txt'
    # outfile_summary = 'Cell_Metadata_summary.txt'
hash_count = 0
norg_count = 0  
nohash_count = 0
nahash_count = 0
duphash_count = 0
hash_read = {}
ecb_hash_dict = {}
c = 0
rg_list = []
hash_list = []
    
header = ['Cell_Barcode', 'ECB_RG', 'RG_color', 'HashTag']
# header = ['Cell_Barcode', 'ECB_RG', 'RG_color', 
#           'HashTag', 'Sample', 'TimePoint', 'Replicate']

with open(os.path.join(out_dir, outfile), 'w') as out,\
     open(os.path.join(out_dir, outfile_min10), 'w') as outmin,\
     open(os.path.join(out_dir, outfile_full), 'w') as outfull,\
     open(os.path.join(out_dir, outfile_summary), 'w')as outsum:
    
    print('\t'.join(header), file=out)
    print('\t'.join(header), file=outmin)
    print('\t'.join(header), file=outfull)
    
    for bc in bc_ref_list:
        rg = 'NA'
        color = 'NA'
        hashtag = 'NA'
        if bc in BC_RG:
            rg = BC_RG[bc]
            color = rg_color[rg]
            if bc in bc_hash:
                hashtag = bc_hash[bc]
                ecb_hash = rg+':'+hashtag
                c += 1
                if hashtag != 'NA' and hashtag != 'Multiplet':
                    outprint = [bc, rg, color, hashtag]
                    print('\t'.join(str(x) for x in outprint), file=out)
                    if rg in RG_cellcount_min10:
                        outminprint = [bc, rg, color, hashtag]
                        print('\t'.join(str(x) for x in outminprint), file=outmin)
                if ecb_hash not in ecb_hash_dict:
                    ecb_hash_dict[ecb_hash] = 1
                else:
                    ecb_hash_dict[ecb_hash] += 1
                if rg not in rg_list:
                    rg_list.append(rg)
                if hashtag not in hash_list:
                    hash_list.append(hashtag)
                
        if bc not in BC_RG:
            norg_count += 1
        
        if bc in bc_hash:
            hashtag = bc_hash[bc]
            if hashtag == 'NA':
                nahash_count += 1
            elif hashtag == 'Multiplet':
                duphash_count += 1
            else:
                hash_count += 1
            
        if bc not in bc_hash:
            nohash_count += 1
            
        if hashtag not in hash_read:
            hash_read[hashtag] = 1
        else:
            hash_read[hashtag] += 1
        
        outfullprint = [bc, rg, color, hashtag]
        print('\t'.join(str(x) for x in outfullprint), file=outfull)
    
    print('Total_BC:', len(bc_ref_list), '\nECB_OK:', 
          len(BC_RG), '\nECB_duplets', multiplet_count,'\nECB_Na', 
          norg_count-multiplet_count, '\nHASH_OK:', hash_count,
          '\nHASH_duplets:', duphash_count, '\nHASH_Na:', nahash_count+nohash_count, 
          '\n\n', file=outsum)
    for k, v in hash_read.items():
        print(k+'\t'+str(v), file=outsum)
    

#### Output ####
# RG_noPrefix.txt 
# Only print if there are >= 30 bc supporting the RG

      
for rg in RG_BC:
    out_name_local = 'RG_' + rg + '.txt'
    if len(RG_BC[rg]) >= 30:
        with open(os.path.join(out_dir, out_name_local), 'w') as outf:
            for bc in RG_BC[rg]:
                print(bc, file=outf)



# outfile_ecbhash = 'Cell_Metadata_ecbs_per_hashSample.txt'

with open(os.path.join(out_dir, outfile_ecbhash), 'w') as outf:

    temprow = []
    print('\t'+'\t'.join(hash_list), file=outf)

    for rg in rg_list:
        temprow = []
        temprow.append(str(rg))
        for hashtag in hash_list:
            tempkey = rg+':'+hashtag
            if tempkey in ecb_hash_dict:
                temprow.append(ecb_hash_dict[str(tempkey)])
            else:
                temprow.append(str(0))
        print('\t'.join([str(x) for x in temprow]), file=outf)

print('Outputs created at %s' %out_dir)
print('Done.')
            


