
from collections import defaultdict
from glob import glob
from os.path import *
from tqdm import tqdm
from subprocess import check_call,check_output

from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import os

# download data
def download(infile,pe=False):
    srrids = open(infile).read().split('\n')
    srrids = [_ for _ in srrids if _]
    for srr in tqdm(srrids):
        if pe:
            cmd = f"fastq-dump --split-files {srr} -O {dirname(infile)}" 
            if not exists(f"{dirname(infile)}/{srr}_1.fastq"):
                check_call(cmd,shell=1)        
        else:
            cmd = f"fastq-dump {srr} -O {dirname(infile)}"  
            if not exists(f"{dirname(infile)}/{srr}.fastq"):
                check_call(cmd,shell=1)        

# get metadata
def get_manifest(indir='./',pe=False,quality=True):
    manifest = f'{realpath(dirname(indir))}/manifest.tab'
    if pe:
        with open(manifest,'w') as f1:
            f1.write('sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n')
            for fq in glob(join(indir,'*_1.fastq')):
                sid = basename(fq).replace('_1.fastq','')
                f1.write(f"{sid}\t{realpath(fq)}\t{realpath(fq).replace('_1.','_2.')}\n")        
    else:
        with open(manifest,'w') as f1:
            f1.write('sample-id\tabsolute-filepath\n')
            for fq in glob(join(indir,'*.fastq')):
                sid = basename(fq).replace('.fastq','')
                f1.write(f'{sid}\t{realpath(fq)}\n')
    return manifest


def get_length_params(in_qzv):
    artifact = Visualization.load(in_qzv)    
    root_dir = '/'.join(str(artifact._archiver.data_dir).split('/')[-2:])
    for_tab_file = join(dirname(in_qzv),root_dir,'forward-seven-number-summaries.tsv')
    rev_tab_file = join(dirname(in_qzv),root_dir,'reverse-seven-number-summaries.tsv')
    cmd = f"cp {in_qzv} {in_qzv.replace('.qzv','zip')}; unzip {in_qzv.replace('.qzv','zip')} -d {dirname(in_qzv)} ;rm -r {in_qzv.replace('.qzv','zip')}"
    if not exists(for_tab_file):
        check_call(cmd,shell=1)

    df = pd.read_csv(for_tab_file,sep='\t',index_col=0)  
    remained_len = df.loc['count',:].max()*0.8
    sub_df = df.loc[:,df.loc["count",:]>remained_len]
    for_r = int(sub_df.columns[-5:].astype(int).values.mean())
    for_l = int(sub_df.columns[0])
    if not exists(rev_tab_file):
        return for_l,for_r
    df = pd.read_csv(rev_tab_file,sep='\t',index_col=0)
    remained_len = df.loc['count',:].max()*0.8
    sub_df = df.loc[:,df.loc["count",:]>remained_len]
    rev_r = int(sub_df.columns[-5:].astype(int).values.mean())
    rev_l = int(sub_df.columns[0])
    
    
    return for_l,for_r,rev_r,rev_l




# srr_id = "ERP012016"
# os.chdir(f'../{srr_id}')


def pe_workflow(srr_id,quick=False,quality=True):
    # os.chdir(f'../{srr_id}')
    indir = f"/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/{srr_id}"
    trimmomatic_dir = "/home-user/thliao/software/Trimmomatic-0.39/"
    if quality:
        for fq in glob(f'{indir}/rawdata/*_1.fastq'):
            sid = basename(fq).replace('_1.fastq','')
            input1 = realpath(fq)
            input2 = realpath(fq).replace('_1.','_2.')
            outdir = dirname(fq).replace('rawdata','QC')
            ofile1 = realpath(input1).replace('/rawdata/','/QC/')
            ofile2 = realpath(input2).replace('/rawdata/','/QC/')
            if not exists(outdir):
                os.makedirs(outdir)
            cmd = f"java -jar /home-user/thliao/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 30 {input1} {input2} {ofile1} {outdir}/{sid}_1.unpaired.fq.gz {ofile2} {outdir}/{sid}_2.unpaired.fq.gz ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"
            if not exists(ofile2) and not exists(ofile1):
                check_call(cmd,shell=1)
            
        manifest = get_manifest(f'{indir}/QC',pe=1)
    else:
        manifest = get_manifest(f'{indir}/rawdata',pe=1)
        
    if quick:
        return
    cmd = f"""qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path {manifest} \
    --output-path {indir}/{srr_id}.qza \
    --input-format PairedEndFastqManifestPhred33V2"""
    if not exists(f"{srr_id}.qza"):
        check_call(cmd,shell=1)


    cmd = f"""qiime demux summarize \
    --i-data {indir}/{srr_id}.qza --o-visualization {indir}/{srr_id}.qzv"""
    if not exists(f"{indir}/{srr_id}.qzv"):
        check_call(cmd,shell=1)



    for_l,for_r,rev_r,rev_l = get_length_params(f"{indir}/{srr_id}.qzv")
    cmd = f"""qiime dada2 denoise-paired \
    --i-demultiplexed-seqs {indir}/{srr_id}.qza \
    --p-trim-left-f {for_l} \
    --p-trim-left-r {rev_l} \
    --p-trunc-len-f {for_r} \
    --p-trunc-len-r {rev_r} \
    --o-representative-sequences {indir}/{srr_id}_repr_dada2.qza \
    --o-table {indir}/{srr_id}_table_dada2.qza \
    --o-denoising-stats {indir}/{srr_id}_stats_dada2.qza \
    --p-n-threads 0 \
    --verbose"""
    if not exists(f"{indir}/{srr_id}_repr_dada2.qza"):
        check_call(cmd,shell=1)


    cmd = f"""
    qiime metadata tabulate \
    --m-input-file {indir}/{srr_id}_stats_dada2.qza \
    --o-visualization {indir}/{srr_id}_stats_dada2.qzv \
    
    qiime feature-table summarize \
    --i-table {indir}/{srr_id}_table_dada2.qza \
    --o-visualization {indir}/{srr_id}_table_dada2.qzv  \
    
    qiime feature-table tabulate-seqs \
    --i-data {indir}/{srr_id}_repr_dada2.qza \
    --o-visualization {indir}/{srr_id}_repr_dada2.qzv
    """
    if not exists(f"{indir}/{srr_id}_repr_dada2.qzv"):
        check_call(cmd,shell=1)


    otutab = f"{indir}/{srr_id}_table_dada2.qza"
    rep = f"{indir}/{srr_id}_repr_dada2.qza"
    stats = f"{indir}/{srr_id}_stats_dada2.qza"
    if not exists(otutab.replace('.qza','.csv')):
        convert2otutab(otutab,
                    otutab.replace('.qza','.csv'))
    if not exists(rep.replace('.qza','.fa')):
        convert2seq(rep,
                rep.replace('.qza','.fa'))
    if not exists(stats.replace('.qza','.csv')):
        convert2otutab(stats,
                    stats.replace('.qza','.csv'))




