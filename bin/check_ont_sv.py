import pysam
file1="/projects/b1042/YueLab/zzhang/ont_pipeline_sample/variants/variants.vcf" #my
file2="/projects/b1042/YueLab/xtwang/workspace/Cancer-nanopore/SCaBER-MD/SCaBER_nanopore.sniffles.vcf" #xiaotao's
def extract_coords_from_sam(sam_file):
    pos_by_chrom={}
    
    with open(sam_file) as fp:
        for line in fp:
            cur_line=line.split('\t')
            if("#" not in line):
                cur_chr,cur_pos=cur_line[0],cur_line[1]
                if(cur_chr) not in pos_by_chrom.keys():
                    pos_by_chrom[cur_chr]=[]
                pos_by_chrom[cur_chr].append(int(cur_pos))
    return pos_by_chrom
my=extract_coords_from_sam(file1)
xt=extract_coords_from_sam(file2)

def calculate_overlap(interval=5):
    total=0
    overlapped=0
    for chromosome, positions in my.items():
        xt_positions=xt[chromosome]
        for cur_pos in positions:
            total+=1
            _,closest_pos=min(enumerate(xt_positions), key=lambda x: abs(x[1]-cur_pos))
            if(abs(closest_pos-cur_pos)<interval+1):
                overlapped+=1
    print("%.1f %%of my SV overlaps with Xiaotao's SV with %dbp interval. "%(((overlapped/total)*100), interval))

calculate_overlap(interval=10)
calculate_overlap(interval=5)
calculate_overlap(interval=0)


