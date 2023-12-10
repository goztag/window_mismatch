#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: goztag
"""

import pysam
import itertools
import pandas as pd
import argparse

#raw vcf file filtered with vcftools  --max-missing 1 --max-alleles 2

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-vcf",dest="vcfinput", type=str, required=True, help="")
    parser.add_argument("-win",dest="windowsize", type=int, required=True, help="")
    parser.add_argument("-total",dest="totalSNPcount", type=int, required=True, help="")
    parser.add_argument("-out",dest="output", type=str, required=True, help="")

args = parser.parse_args()

vcfile=args.vcfinput
out=args.output


window_size = int(args.windowsize)  # define window size
total=int(args.totalSNPcount) # define total number of SNPs


vcf = pysam.VariantFile(vcfile)


samples = list(vcf.header.samples)

# pairs combinations of individuals
individual_pairs = list(itertools.combinations(samples, 2))


# function to calculate number of mismatches per pair
def calculate_mismatches(record, ind1, ind2):
    
    gt1=record.samples[ind1]['GT']
    gt2=record.samples[ind2]['GT']
    
    gtsum1=sum(gt1)
    gtsum2=sum(gt2)

    mismatch=abs(gtsum1 - gtsum2)
    
    return mismatch



def meanmismatch_perwindow():
    
    # open a dictionary for the mismatch values per window r:pair & c:meanmismatch

    pairwise_mismatches = {pair: [] for pair in individual_pairs}


    current_window_start = 0
    current_window_end = current_window_start + window_size

    #define window size for loop use separately, 
    #as it may change according to the final window size
    window_size2=int(args.windowsize)

    count=1 #snp count
    
    #dictionary to temporarily store current window mismatches per pair
    pair_mm={pair: [] for pair in individual_pairs} 


    for record in vcf.fetch(): #read vcf lines
        
        if count <= current_window_end : #if the snp is within the current window
        
            for pair in individual_pairs: #calculate mismatches for that snp
                
                ind1=pair[0]
                ind2=pair[1]
            
                mm=calculate_mismatches(record, ind1, ind2)

                pair_mm[pair].append(mm)
        
           
            
        if count > current_window_end:  #if the snp is outside the current window
            
           
            for pair in individual_pairs:
                
                ind1=pair[0]
                ind2=pair[1]
            
                #calculate mean mismatch for current window
                mean_mismatch = sum(pair_mm[pair]) / window_size2

                pairwise_mismatches[pair].append(mean_mismatch)
            
            #next window 
            current_window_start = current_window_end
            
            
            #if not in the last window
            if current_window_end + window_size <= total :
            
                current_window_end += window_size
                
            elif current_window_end + window_size > total : #if in last window
                
                #find the size of last window 
                window_size2= total - current_window_end 
                #assign as last window
                current_window_end = total
                
            #reset for next window    
            pair_mm={pair: [] for pair in individual_pairs}

            #current snp is the 1st snp of next window, so calculate the mismatch
            #and store, this step is only for indices of window ends
            for pair in individual_pairs:
                
                ind1=pair[0]
                ind2=pair[1]
            
                mm=calculate_mismatches(record, ind1, ind2)

                pair_mm[pair].append(mm)
        #next snp    
        count += 1 

    #only for the last window
    #since the previous loop ends prior to the mean calculation

    for pair in individual_pairs:
        
        ind1=pair[0]
        ind2=pair[1]

        mean_mismatch = sum(pair_mm[pair]) / window_size2

        pairwise_mismatches[pair].append(mean_mismatch)
        
    
    return pairwise_mismatches

    
mismatch_result=meanmismatch_perwindow()


# Convert to data frame
pairwise_mismatches_df = pd.DataFrame.from_dict(mismatch_result, orient='index')

# if desired,rows: windows and columns: pairs
pairwise_mismatches_df = pairwise_mismatches_df.transpose()

#csv
pairwise_mismatches_df.to_csv(out, index=True)






