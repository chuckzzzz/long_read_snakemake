# egrep -v "^#" variants.vcf > variants_no_header.vcf
# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
import io
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


plt.rcParams["figure.figsize"] = (20,15)


input_file="/projects/b1042/YueLab/zzhang/ont_pipeline_sample/variants/variants.vcf"
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
df=read_vcf(input_file)

# extract number of reads supporting a variant
info=df["INFO"].to_numpy()
RE_list=[int(l.split(";")[RE_index].strip("RE=")) for l in info]
# calculate distribution of variant types
SV_list=[l.split(";")[SV_TYPE_index].split("SVTYPE=")[1] for l in info]
# get SV position info
chr_list=df["CHROM"]
pos_list=df["POS"]

analysis_df=pd.DataFrame({
    "CHROM":chr_list,
    "POS":pos_list,
    "RE":RE_list,
    "SV":SV_list
})

fig, axes = plt.subplots(2)
sns.histplot(data=analysis_df["SV"],ax=axes[0])
sns.histplot(data=analysis_df["RE"],ax=axes[1],bins=30)
plt.savefig("sv_vis.png")

# calculate 99 percentile RE
thres = np.percentile(analysis_df["RE"], 99) # return 50th percentile, e.g median.
target_df=analysis_df[analysis_df['RE']>thres]
target_df=target_df.reset_index()
target_df.to_csv("targets.tsv",sep="\t",index=False)