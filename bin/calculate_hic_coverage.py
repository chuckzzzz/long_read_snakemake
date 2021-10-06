'''
Calculate the coverage of HiC contact map. The script is an adaptation on https://github.com/XiaoTaoWang/NeoLoopFinder/blob/master/neoloop/cnv/runcnv.py

Only works with mcool file at the moment. For single cool file, cal get_marginals directly by passing cool file patt. 
'''

import os, cooler, pyBigWig
from cooler import balance
import numpy as np
import pandas as pd
from collections import Counter
import argparse
import h5py
import pathlib

def get_marginals(uri, exclude=['M', 'Y', 'MT', 'EBV'], chunksize=int(1e7), nproc=1):

    clr = cooler.Cooler(uri)

    if nproc > 1:
        pool = balance.Pool(nproc)
        map_ = pool.imap_unordered
    else:
        map_ = map

    nnz = clr.info['nnz']
    n_bins = clr.info['nbins']
    edges = np.arange(0, nnz+chunksize, chunksize)
    spans = list(zip(edges[:-1], edges[1:]))

    marg = (
        balance.split(clr, spans=spans, map=map_, use_lock=False)
            .prepare(balance._init)
            .pipe([])
            .pipe(balance._marginalize)
            .reduce(balance.add, np.zeros(n_bins))
    )
    table = clr.bins()[:][['chrom', 'start', 'end']]
    table['Coverage'] = marg.astype(int)
    pool = []
    chroms = [c for c in clr.chromnames if ((not c.lstrip('chr') in exclude) and (not '_' in c))]
    for chrom in chroms:
        pool.append(table[table['chrom']==chrom])
    
    table = pd.concat(pool)

    return table

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Calculate Hi-C coverage. ')
    parser.add_argument('-c', help='Path to hi-c file',required=True)
    parser.add_argument('-o', help='Path to output directory',required=True)

    args = parser.parse_args()
    # path to mcool file
    hic_file=args.c
    # path to output directory
    output_prefix=args.o 
    pathlib.Path(output_prefix).mkdir(parents=True, exist_ok=True) 

    # unpack cool files from mcool file
    h5=h5py.File(hic_file,'r')
    all_resolutions=list(h5['resolutions'].keys())
    cool_files=[hic_file+"::/resolutions/"+reso for reso in all_resolutions]

    # apply the coverage calculation function to each file sequentially
    coverage_tables=list(map(get_marginals,cool_files))
    for table in coverage_tables:
        cur_bin=str(table.loc[0,"end"]-table.loc[0,"start"])
        output_path=os.path.join(output_prefix,"reso_"+cur_bin+".tsv")
        table.to_csv(output_path,sep='\t',index=False)
# python calculate_hic_coverage.py -c /projects/b1042/YueLab/xtwang/ziyang/SCABER-Arima-allReps-filtered.hg38.mcool -o /projects/b1042/YueLab/zzhang/ont_pipeline_sample/analysis/hic_coverage
