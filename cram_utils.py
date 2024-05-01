import pysam
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import edlib
import csv
import json
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import re
from multiprocessing import Pool
import mappy

UP_seq = "TCTTCAGCGTTCCCGAGA"
UP_seq_revcomp = "TCTCGGGAACGCTGAAGA"

const_3prime_full_forward = 'TCAGACGTGTGCTCTTCCGATCT' #right
const_3prime_full_reverse = 'AGATCGGAAGAGCACACGTCTGA' #right
const_5prime_full_forward = 'CTACACGACGCTCTTCCGATCT' # left

right_const = const_3prime_full_reverse[:20]
left_const = const_5prime_full_forward[-20:]

N_read_extract = 1000000

print(N_read_extract)


def cram_fastq_split_by_lines(indir, sample, cores, lines=4e6):
    cram_file = f"{indir}/{sample}.cram"
    crai_file = f"{indir}/{sample}.cram.crai"

    if os.path.isfile(crai_file):
        print(f"cram {crai_file} index exists")
    else:
        print(f"cram {crai_file} index does not exit, will index")

        cram_index_command = f"samtools index -@{int(cores)} {cram_file}"

        subprocess.call(cram_index_command, shell=True)

    splitted_file = f"{indir}/{sample}/split/{sample}.part_000.fastq"

    if os.path.isfile(splitted_file):
        print(splitted_file, " splitted fastq exists, skip splitting")
    else:
        print(splitted_file, " splitted fastq does not exist")

        split_dir = f"{indir}/{sample}/split"
        if not os.path.exists(split_dir):
            os.makedirs(split_dir)
            print(f"{split_dir} created")
        else:
            print(f"{split_dir} already exists")

        split_name_pattern = f"{split_dir}/{sample}.part_"

        # convert to fastq with samtools and split by n=lines (4x number of reads)
        # add a suffix with 3 digits and prefix of 'split_name_pattern'
        cram_to_fastq_command = (
            f"samtools fastq -@{int(cores)} {cram_file} | "
            f"split -d -a 3 -l {int(lines)} "
            f"--additional-suffix=.fastq -  {split_name_pattern}"
        )

        subprocess.call(cram_to_fastq_command, shell=True)


def seq_counter(seq_dict, seq_instance):
    if seq_dict.get(seq_instance) is None:
        seq_dict[seq_instance] = 1
    else:
        seq_dict[seq_instance] += 1


def quad_dict_store(quad_dict, quad_key, quad_items):
    if quad_dict.get(quad_key) is None:
        quad_dict[quad_key] = [quad_items]
    else:
        quad_dict[quad_key].extend([quad_items])


def edit_match(input_seq, target_seq, max_dist):
    if input_seq == target_seq:
        dist = 0
        match = True
    else:
        edit = edlib.align(input_seq, target_seq, "NW", "path", max_dist)
        dist = edit["editDistance"]
        if dist >= 0 and dist <= max_dist:
            cigar = edit["cigar"]
            if "D" in cigar or "I" in cigar:
                match = False
                dist = "indel"
            else:
                match = True
        else:
            match = False

    return (match, dist)


def quality_calc(seq, quals, bases_dict, quals_dict):
    for i in range(len(seq)):
        if bases_dict.get(str(i)) is None:
            bases_dict[str(i)] = {}
            seq_counter(bases_dict[str(i)], seq[i])
        else:
            seq_counter(bases_dict[str(i)], seq[i])

        if quals_dict.get(str(i)) is None:
            quals_dict[str(i)] = {}
            seq_counter(quals_dict[str(i)], quals[i])
        else:
            seq_counter(quals_dict[str(i)], quals[i])


def quality_df(quals_dict):
    quals_df = pd.DataFrame(quals_dict)
    quals_df = quals_df.T
    quals_df = quals_df.fillna(0)
    quals_df = quals_df.stack()
    quals_df = quals_df.reset_index()
    # quals_df.columns = ['base', 'quality', 'tot_count']
    # quals_df['mult'] = quals_df.quality * quals_df.tot_count
    # quals_df_grouped = quals_df.groupby('base').sum()
    quals_df.columns = ["position", "base_qual", "tot_count"]
    quals_df.position = quals_df.position.astype("int")
    # quals_df[quals_df.position.isin(np.arange(10))]
    counts_df = quals_df.groupby("position").sum()
    quals_df["position_cnt"] = quals_df.position.apply(
        lambda x: counts_df.loc[x].tot_count
    )
    quals_df["freq"] = quals_df.tot_count / quals_df.position_cnt * 100
    return quals_df


def extract_trimmed_fastq_pairs(indir, sample, part, limit):
    i = 0
    max_dist = 3

    ultima_fastq = f"{indir}/{sample}/split/{sample}.part_{part}.fastq"

    R1_fastq = f"{indir}/{sample}/split/{sample}_R1.part_{part}.fastq"
    R2_fastq = f"{indir}/{sample}/split/{sample}_R2.part_{part}.fastq"

    r_qual_dict = {}
    r_base_dict = {}

    do_qc = True

    #if os.path.isfile(bcs_json):
    #    print(bcs_json, " exists, skip")
    #    return

    R1 = open(R1_fastq, "w")
    R2 = open(R2_fastq, "w")

    with pysam.FastxFile(ultima_fastq) as R:
        for r in tqdm(R):
            
            i += 1

            seq = r.sequence

            rlen = len(seq)
            
            if rlen>130 and rlen<170:
            
                end_seq = seq[-50:]
                accept_r2 = False
                pos_con_in_end = end_seq.find(right_const)
                if pos_con_in_end>=0:
                    accept_r2 = True
                    dist = 0
                else:
                    edit = edlib.align(right_const, end_seq,'HW','locations',3)
                    dist = edit['editDistance']
                    if dist>=0:
                        accept_r2 = True
                        locs = edit['locations'][0]
                        pos_con_in_end = locs[0]

                begin_seq = seq[:50]
                accept_r1 = False
                pos_con_in_begin = begin_seq.find(left_const)
                if pos_con_in_begin >=0 :
                    accept_r1 = True
                    dist_left_const = 0
                    pos_con_in_begin += len(left_const)
                else:
                    edit = edlib.align(left_const, begin_seq, 'HW', 'locations', 3)
                    dist = edit['editDistance']
                    if dist >=0 :
                        accept_r1 = True
                        locs = edit['locations'][0]
                        pos_con_in_begin = locs[1]+1
                        dist_left_const = dist

                if accept_r2 and accept_r1:
                    
                    qual = r.quality
                    trim_begin  = pos_con_in_begin
                    r1_seq = seq[trim_begin:trim_begin+50]
                    r1_qual = qual[trim_begin:trim_begin+50]
                    
                    trim_end  = 50 - pos_con_in_end
                    r2_seq = mappy.revcomp(seq[-trim_end-50:-trim_end])
                    r2_qual = qual[-trim_end-50:-trim_end][::-1]
                    
                    R1.write(f"@{r.name}_1\n")
                    R1.write(f"{r1_seq}\n")
                    R1.write("+\n")
                    R1.write(f"{r1_qual}\n")

                    R2.write(f"@{r.name}_2\n")
                    R2.write(f"{r2_seq}\n")
                    R2.write("+\n")
                    R2.write(f"{r2_qual}\n")
                    
                    #if do_qc and i % 500 == 0:
                    #    
                    #    quals = r.get_quality_array()
                    #    quality_calc(seq, quals, r_base_dict, r_qual_dict)
                        
                    
                    
                    #r2_all.append(r2_seq)
                    #print('r2',pos_con_in_end, mappy.revcomp(seq[-trim_end-50:-trim_end]))

                    
                    #r1_all.append(r1_seq)
                    #print('r1',pos_con_in_begin, r1_seq,dist_left_const)
                
            
            if i > N_read_extract and limit:
                break

    #r_qual_df = quality_df(r_qual_dict)
    #r_base_df = quality_df(r_base_dict)
    #r_qual_df.to_csv(ultima_fastq.replace(".fastq", "_quals.csv"))
    #r_base_df.to_csv(ultima_fastq.replace(".fastq", "_bases.csv"))

    R1.close()
    R2.close()

                
                
def seq_slice(read_seq, bc_intervals, umi_intervals):
    bc = "".join([read_seq[intv[0] : intv[1]] for intv in bc_intervals])
    umi = "".join([read_seq[intv[0] : intv[1]] for intv in umi_intervals])

    return (bc, umi)


def find_sub_fastq_parts(indir, sample):
    #pattern = re.compile(r"_R1.part_(.*?)\.fastq")
    pattern = re.compile(r"^(?!.*R[12]).*part_\d{3}\.fastq")
    print(pattern)
    all_files = os.listdir(f"{indir}/{sample}/split/")
    # part + 3 digits because we did split suffix with 3 digits
    all_parts = [f.split(".part_")[1][:3] for f in all_files if pattern.search(f)]
    parts = sorted(np.unique(all_parts))

    return parts