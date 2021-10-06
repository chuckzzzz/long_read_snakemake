'''
Plot Hi-C plots at different resolutions for candidate SV identified by Sniffles. 
'''

import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import h5py
import cooler
import pathlib
import os
import argparse
import sys

# plot parameters for a 3X3 grid plot
plt.rcParams["figure.figsize"] = (24, 24)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze SV results and generate visualizations.')
    parser.add_argument('-s', help='Path to strunctural variants file', required=True)
    parser.add_argument('-c', help='Path to hi-c file',required=True)
    parser.add_argument('-o', help='Path to output directory. Default is current directory. ', default="./")
    parser.add_argument('-b', help='Total bins around SV to plot. Default is 400. ', default=400)
    parser.add_argument('-m', help='Flag for mcool format input. Default is false. ', dest='m', action='store_true',default=False)
    parser.add_argument('-f', help='Flag for scale. If true will remove columns and rows in hi-c matrix that only has nan(no signal). Default is false. ', dest='f', action='store_true',default=False)
    parser.add_argument('-n', action='store_true', default=False, dest='n', help='Flag for whether to normalize matrix. Default is false. ')

    args = parser.parse_args()
    # flag for whether input file is an mcool file
    mcool=args.m
    # number of bins to show in the visualization
    num_bins=args.b 
    # target SV path
    target_sv_path=args.s 
    # hic cool format file path
    hic_file=args.c
    # output directory, will be created if doesn't exist
    output_prefix=args.o 
    pathlib.Path(output_prefix).mkdir(parents=True, exist_ok=True) 
    # scale or not
    scale=args.f
    # normalize or not
    balance=args.n
    # check if required inputs are valid
    assert os.path.exists(target_sv_path), "SV target file could not be located!"
    assert os.path.exists(hic_file), "hi-C target file could not be located!"

    # parse target SV file
    sv_df=pd.read_csv(target_sv_path,sep='\t')

    if(mcool):
        
        # retreive all resolutions in the mcool file
        h5=h5py.File(hic_file,'r')
        all_resolutions=list(h5['resolutions'].keys())
        cool_files=[hic_file+"::/resolutions/"+reso for reso in all_resolutions]
        for index, row in sv_df.iterrows():
            # iterate through the candidates
            cur_chr,cur_pos,cur_RE,cur_SV=row['CHROM'],int(row['POS']),row['RE'],row['SV']

            fig, axes = plt.subplots(nrows=3, ncols=3)

            # iterate through different resolutions
            for idx,cur_cool in enumerate(cool_files):
                c=cooler.Cooler(cur_cool)
                bin_size=c.info["bin-size"]
                chrom_sizes=pd.DataFrame(c.chromsizes)
                cur_chr_max=int(chrom_sizes.loc[cur_chr])

                rounding_factor=-1*(int(np.log10(bin_size)))

                # round by bin size
                cur_pos_rounded=round(cur_pos,rounding_factor)
                
                # calculate start and end coordinates for the Hi-C map
                cur_start,cur_end=cur_pos_rounded-(num_bins//2)*bin_size, cur_pos_rounded+(num_bins//2)*bin_size
                
                # handle out of boundary cases
                cur_start=max(1,cur_start)
                cur_end=min(cur_end,cur_chr_max)
    
                # retreive plotting matrix 
                fetching_string=cur_chr+":"+str(cur_start)+"-"+str(cur_end)
                target_mtx=c.matrix(balance=balance).fetch(fetching_string)

                # debug
                # t1=c.matrix(balance=True).fetch(fetching_string)
                # t2=c.matrix(balance=False).fetch(fetching_string)
                # print(np.allclose(t1, t2, equal_nan=True))

                # prevent nans in the matrix after normalization
                target_mtx=np.log10(target_mtx)
                if(scale):
                    target_mtx=target_mtx[:,~np.all(np.isinf(target_mtx), axis=0)]
                    target_mtx=target_mtx[:,~np.all(np.isnan(target_mtx), axis=0)]
                    target_mtx=target_mtx[~np.all(np.isinf(target_mtx), axis=1),:]
                    target_mtx=target_mtx[~np.all(np.isnan(target_mtx), axis=1),:]
                
                # calculate subplot index and plot 
                row_idx=(idx+1)%3 # can be adjusted later
                col_idx=(idx+1)//3
                im=axes[row_idx-1,col_idx-1].matshow(target_mtx, cmap='YlOrRd')
                plt.colorbar(im,ax=axes[row_idx-1,col_idx-1])
                axes[row_idx-1,col_idx-1].title.set_text(fetching_string+" | bin: %d"%bin_size)

            sv_id="_".join([cur_SV,cur_chr,str(cur_pos)])
            if(balance):
                output_path=os.path.join(output_prefix,sv_id+".normalized.png")
            else:
                output_path=os.path.join(output_prefix,sv_id+".png")
            cur_main_title="SV Type: %s | Supporting Reads: %s" %(cur_SV,cur_RE)
            plt.suptitle(cur_main_title,fontsize=28)
            plt.savefig(output_path)
            plt.clf()
            plt.close()
    else:
        sys.exit("Single file function not implemented yet. ")    
        


# python plot_hic_graph.py -s /projects/b1042/YueLab/zzhang/ont_pipeline_sample/analysis/target_sv.tsv -c /projects/b1042/YueLab/xtwang/ziyang/SCABER-Arima-allReps-filtered.hg38.mcool -o /projects/b1042/YueLab/zzhang/ont_pipeline_sample/analysis/hic -m -f

# python plot_hic_graph.py -s /projects/b1042/YueLab/zzhang/ont_pipeline_sample/analysis/sv_to_zoom_in.tsv -c /projects/b1042/YueLab/xtwang/ziyang/SCABER-Arima-allReps-filtered.hg38.mcool -o /projects/b1042/YueLab/zzhang/ont_pipeline_sample/analysis/hic_zoom_in -m

# python plot_hic_graph.py -s /projects/b1042/YueLab/zzhang/ont_pipeline_sample/analysis/sv_to_zoom_in.tsv -c /projects/b1042/YueLab/xtwang/ziyang/SCABER-Arima-allReps-filtered.hg38.mcool -o /projects/b1042/YueLab/zzhang/ont_pipeline_sample/analysis/hic_zoom_in -m -b