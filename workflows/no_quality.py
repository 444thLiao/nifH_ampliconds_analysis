

from subprocess import check_call
from os.path import *
from Bio import SeqIO
from glob import glob
from .q2function import convert2otutab,convert2seq
from tqdm import tqdm

def no_q_workflow(sra):
    indir = f"/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/{sra}"
    
    all_records = []
    for fq in glob(join(indir,'rawdata','*.fastq')):
        srr = basename(fq).replace('.fastq','')
        records = SeqIO.parse(fq,'fastq')
        for r in records:
            r.name = r.description = ''
            r.id = f"{srr}_{r.id}"
            all_records.append(r)
    with open(join(indir,'all.fna'),'w') as f1:
        SeqIO.write(all_records,f1,'fasta-2line')
        
        
    cmd = f"""qiime tools import \
  --input-path {indir}/all.fna \
  --output-path {indir}/all.qza \
  --type 'SampleData[Sequences]'
  """   
    if not exists(f"{indir}/all.qza"):
        check_call(cmd,shell=1)
        
    cmd = f"""qiime vsearch dereplicate-sequences \
  --i-sequences {indir}/all.qza \
  --o-dereplicated-table {indir}/{sra}_dep_table.qza \
  --o-dereplicated-sequences {indir}/{sra}_dep.qza \
    --verbose
    """
    if not exists(f"{indir}/{sra}_dep.qza"):
        check_call(cmd,shell=1)
        
        
    cmd = f"""qiime vsearch cluster-features-de-novo \
  --i-table {indir}/{sra}_dep_table.qza \
  --i-sequences {indir}/{sra}_dep.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table {indir}/{sra}_table_dada2.qza \
  --o-clustered-sequences {indir}/{sra}_repr_dada2.qza \
  --verbose
    """
    if not exists(f"{indir}/{sra}_repr_dada2.qza"):
        check_call(cmd,shell=1)
        
        
    
    otutab = f"{indir}/{sra}_table_dada2.qza"
    rep = f"{indir}/{sra}_repr_dada2.qza"
    stats = f"{indir}/{sra}_stats_dada2.qza"
    if not exists(otutab.replace('.qza','.csv')):
        convert2otutab(otutab,
                    otutab.replace('.qza','.csv'))
    if not exists(rep.replace('.qza','.fa')):
        convert2seq(rep,
                rep.replace('.qza','.fa'))
    
    vsearch = "/home-user/thliao/anaconda3/envs/qiime2-2020.8/bin/vsearch"
    infa = 'all.fna'
    filtered_fa = f'{indir}/filtered.fna'
    derep_fa = f'{indir}/derep.fna'
    derep_uc = f'{indir}/derep.uc'
    sorted_d2_fa = f'{indir}/derep_d2.fna'
    cluster_uc = f'{indir}/cluster.uc'
    rep_fa = f'{indir}/{sra}_repr_dada2.fa'
    map_output = f'{indir}/final_map.uc'
    raw_otutab = f'{indir}/{sra}_table_dada2.csv'
    
    
    # denovo_nonchimer_fa = f'{indir}/denovo_nonchimer.fna'
    
    # cmd = f"{vsearch} --fastx_filter {indir}/{infa} --fastq_maxee 1.0 --fastq_trunclen 240 --fastaout {filtered_fa}"
    # if not exists(filtered_fa):
    #     check_call(cmd,shell=1)
    cmd = f"{vsearch} --derep_fulllength {indir}/{infa} --output {derep_fa} -sizeout --fasta_width 0 --uc {derep_uc}"
    if not exists(derep_fa):
        check_call(cmd,shell=1)
    cmd = f"{vsearch} --sortbysize  {derep_fa} --output {sorted_d2_fa} --minsize 2"
    if not exists(sorted_d2_fa):
        check_call(cmd,shell=1)

    cmd = f"{vsearch} --cluster_size {sorted_d2_fa} --id 0.97 --sizein --sizeout --fasta_width 0 --uc {cluster_uc} --relabel OTU --centroids {rep_fa}"
    if not exists(rep_fa):
        check_call(cmd,shell=1)

    cmd = f"{vsearch} --usearch_global {indir}/{infa} --db {rep_fa} --strand plus --id 0.97 --uc {map_output} --otutabout {raw_otutab}"
    if not exists(raw_otutab):
        check_call(cmd,shell=1)

    





import gzip
all_records = []

with open(join('./','all.fna'),'w') as f1:
    for fq in tqdm(glob('./f8c9eb7f-ff53-4737-96ba-de172055a16c/data/*.fastq.gz')):
        srr = basename(fq).split('_')[0]
        records = SeqIO.parse(gzip.open(fq,'rt'),'fastq')
        for r in records:
            r.name = r.description = ''
            r.id = f"{srr}_{r.id}"
            # all_records.append(r)
            SeqIO.write(r,f1,'fasta-2line')