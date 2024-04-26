import argparse
from multiprocessing import Pool
import cram_utils
import os

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cores', type=str)
parser.add_argument('-i', '--indir', type=str)
parser.add_argument('-s', '--sample', type=str)
parser.add_argument('-l', '--limit', default=False, action='store_true')

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
limit = args.limit

######################################################

cram_utils.cram_fastq_split_by_lines(indir, sample, cores, 4e7)

######################################################

"""
parts = bc_umi_utils.find_sub_fastq_parts(indir,sample)
args = [(indir, sample, part, limit, read1_struct, read2_struct) for part in parts]

pool = Pool(int(cores))
results = pool.starmap(bc_umi_utils.extract_bc_umi_dict, args)
pool.close()
pool.join()

######################################################

bc_umi_utils.aggregate_dicts(indir,sample,'anchors')
bc_umi_utils.aggregate_dicts(indir,sample,'targets')

bc_umi_utils.aggregate_stat_dicts(indir,sample,'adapter_edits')

######################################################

qc_pdf_file = f'{indir}/{sample}/{sample}_QC.pdf'

if os.path.isfile(qc_pdf_file):
    print(qc_pdf_file,' exists, skip')
else:
    qc_pdfs = PdfPages(qc_pdf_file)
    bc_umi_utils.whitelist_rankplot(indir,sample,'anchors',qc_pdfs,max_anchors)
    bc_umi_utils.whitelist_rankplot(indir,sample,'targets',qc_pdfs,max_targets)
    qc_pdfs.close()
    

"""