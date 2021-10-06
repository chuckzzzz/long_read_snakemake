'''
Query a coverage table to with a ucsc format string. Stores the resulting table at designated location.  
'''

import pandas as pd
import argparse
import pathlib
import os

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Query a coverage table')
    parser.add_argument('-c', help='Path to coverage file',required=True)
    parser.add_argument('-s', help='Query String in UCSC genome browser format. Example would be chr1:142213353-143413353',required=True)
    parser.add_argument('-o', help='Path to output directory',required=True)

    args = parser.parse_args()
    # path to mcool file
    coverage_file=args.c
    # path to output directory
    output_prefix=args.o 
    pathlib.Path(output_prefix).mkdir(parents=True, exist_ok=True) 
    # parse the query string
    query_string=args.s 
    cur_chr=query_string.split(':')[0]
    start,end=query_string.split(':')[1].split('-')
    start,end=int(start),int(end)
    
    # subset dataframe by chrosomome then by position
    all_df=pd.read_csv(coverage_file,sep='\t')
    all_df=all_df.loc[lambda all_df:all_df['chrom']==cur_chr]
    output_df=all_df[(all_df['start']>=start) & (all_df['end']<=end)]
    
    # store output
    output_file=os.path.join(output_prefix,query_string+".tsv")
    output_df.to_csv(output_file,index=False,sep='\t')
