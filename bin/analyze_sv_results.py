# egrep -v "^#" variants.vcf > variants_no_header.vcf
# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
import io
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import pathlib


plt.rcParams["figure.figsize"] = (20,15)
RE_index=-4
SV_TYPE_index=-9

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

#input_file="/projects/b1042/YueLab/zzhang/ont_pipeline_sample/variants/variants.vcf"
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze SV results and generate visualizations.')
    parser.add_argument('-i', help='Path to input file', required=True)
    parser.add_argument('-o', help='Path to output directory', default="./")
    parser.add_argument('-p', help='Percentile of supporting reads number', default=99)

    args = parser.parse_args()

    # parse the input arguments
    input_file=args.i
    re_percentile=args.p
    output_prefix=args.o
    output_sv_tsv=os.path.join(output_prefix,"target_sv.tsv")
    output_vis_fig=os.path.join(output_prefix,"sv_vis.png")

    # creates the output directory recursively if it does not exsit
    pathlib.Path(output_prefix).mkdir(parents=True, exist_ok=True) 
    
    # read vcf file
    df=read_vcf(input_file)

    # extract number of reads supporting a variant
    info=df["INFO"].to_numpy()
    RE_list=[int(l.split(";")[RE_index].strip("RE=")) for l in info]
    # calculate distribution of variant types
    SV_list=[l.split(";")[SV_TYPE_index].split("SVTYPE=")[1] for l in info]
    # get SV position info
    chr_list=df["CHROM"]
    pos_list=df["POS"]

    # construct a separate dataframe for easy visualization 
    analysis_df=pd.DataFrame({
        "CHROM":chr_list,
        "POS":pos_list,
        "RE":RE_list,
        "SV":SV_list
    })

    fig, axes = plt.subplots(2)
    sns.histplot(data=analysis_df["SV"],ax=axes[0])
    sns.histplot(data=analysis_df["RE"],ax=axes[1],bins=30)
    plt.savefig(output_vis_fig)

    # calculate 99 percentile RE
    thres = np.percentile(analysis_df["RE"], re_percentile) 
    target_df=analysis_df[analysis_df['RE']>thres]
    target_df=target_df.reset_index(drop=True)
    target_df.to_csv(output_sv_tsv,sep="\t",index=False)


# python analyze_sv_result
