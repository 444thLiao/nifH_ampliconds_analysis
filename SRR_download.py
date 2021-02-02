"""
Download those srr
"""
from os.path import *
import pandas as pd
from subprocess import check_call
from tqdm import tqdm
import os
import warnings
warnings.filterwarnings("ignore")

os.chdir("./amplicons_study")

header = ['qaccver',
          'saccver',
          'pident',
          'length',
          'mismatch',
          'gapopen',
          'qstart',
          'qend',
          'sstart',
          'send',
          'evalue',
          'bitscore']


def check_done(infile, pe=False):
    # check_download
    srrids = open(infile).read().split('\n')
    srrids = [_ for _ in srrids if _]
    if not srrids:
        return False
    sra_id = infile.split('/')[-3]

    wrong_layout = False
    for srr in srrids:
        if pe:
            # if pe have not complete srr_download,return False
            if not exists(f"{dirname(infile)}/{srr}_1.fastq"):
                return False
            # if set pe, but can't detect _2.fastq
            # this must be wrong layout set.
            elif exists(f"{dirname(infile)}/{srr}_1.fastq") and (not exists(f"{dirname(infile)}/{srr}_2.fastq")):
                os.system(
                    f"mv {dirname(infile)}/{srr}_1.fastq {dirname(infile)}/{srr}.fastq")
                wrong_layout = True
            # if set pe, but only detect .fastq
            elif exists(f"{dirname(infile)}/{srr}.fastq"):
                wrong_layout = True
                break
        else:
            # if se set
            if (not exists(f"{dirname(infile)}/{srr}.fastq")) and (not exists(f"{dirname(infile)}/{srr}_1.fastq")):
                return False
            elif (exists(f"{dirname(infile)}/{srr}_1.fastq")) and (exists(f"{dirname(infile)}/{srr}_2.fastq")):
                wrong_layout = True
                print(f"{sra_id} should be paried")
                break
    if wrong_layout:
        print(f"{sra_id} should be Single")
        return False
    # if complete qiime2, return False
    if exists(f"/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/{sra_id}/{sra_id}_repr_dada2.fa"):
        return False
    return True


def download(infile, pe=False, check=False):
    srrids = open(infile).read().split('\n')
    srrids = [_ for _ in srrids if _]
    sra_id = infile.split('/')[-3]
    # f"fastq-dump --split-files {srr} -O {dirname(infile)}"
    for srr in tqdm(srrids):
        cmd = f"fastq-dump {srr} -O {dirname(infile)} "
        cmd += ' --split-files' if pe else ''
        if pe and exists(f"{dirname(infile)}/{srr}_1.fastq"):
            continue
        elif (not pe) and exists(f"{dirname(infile)}/{srr}.fastq"):
            continue
        else:
            pass
        if not check:
            try:
                check_call(cmd, shell=1)
            except:
                pass
        else:
            print(f'{srr} of {sra_id} need download')
            return sra_id


df2 = pd.read_excel(
    './data/SraRunTable.xlsx')
df = pd.read_excel(
    './data/nifH metagenome_filtered_16S.xlsx')

# df2.loc[df2['Experiment'].isin(set(df['Experiment Accession']))]

removed_srr = {'SRP135966': ['SRR6854367'],
               "ERP108575": ["ERR2560196"]}


gb = df2.groupby('SRA Study')
all_list = {}
for study_acc, _gidx in gb.groups.items():
    sub_df = df2.loc[_gidx, :]
    sub_df = sub_df.loc[~sub_df['AvgSpotLen'].isna(), :]
    srr_ids = sub_df['Run']

    if not exists(f'./{study_acc}/rawdata'):
        os.system(f'mkdir -p ./{study_acc}/rawdata')

    if study_acc in removed_srr:
        srr_ids = [_ for _ in srr_ids if _ not in removed_srr[study_acc]]
    with open(f'./{study_acc}/rawdata/srr_ids.list', 'w') as f1:
        f1.write('\n'.join(srr_ids))
    # all_list[f'{study_acc}/rawdata/srr_ids.list'] = list(sub_df['LibraryLayout'])[0]
    all_list[study_acc] = list(sub_df['LibraryLayout'])[0]

run2info = df2.set_index('Run')
run2info = run2info.to_dict(orient='index')

removed_sra = ["ERP020550",  # not nifH
               "ERP108575",  # need barcode to demux
               "ERP015983",  # need barcode to demux
               "SRP095888",  # not nifH 16S
               #"SRP053025", # tmp, (only some of them are nifH)
               "SRP139936",  # no usable reads for processed

               "SRP095769",  # error during dada2 denoise , reverse pair show weird quality score, maybe try it without the reverse parts
               ]
no_quality_ids = [
    "SRP227487",  # no quality information
    "SRP273539",   # no quality information
    "SRP273532",   # no quality information
    "SRP089934",   # no quality information
    "SRP215034"   # no quality information
]
removed_sra += no_quality_ids
all_list['SRP053025'] = 'SINGLE'  # liar !!!!
all_list['SRP106861'] = 'SINGLE'  # liar !!!!
all_list['SRP108172'] = 'SINGLE'  # liar !!!!
all_list['SRP125896'] = 'SINGLE'  # liar !!!!
all_list['SRP128014'] = 'SINGLE'  # liar !!!!
all_list['SRP153053'] = 'SINGLE'  # liar !!!!
all_list['SRP158049'] = 'SINGLE'  # liar !!!!
all_list['SRP051798'] = 'SINGLE'  # liar !!!!
all_list['SRP113682'] = 'SINGLE'  # liar !!!!
all_list["SRP078449"] = 'SINGLE'
all_list["SRP171566"] = 'SINGLE'
all_list["SRP223660"] = 'SINGLE'
all_list["SRP273532"] = 'SINGLE'
all_list["SRP273539"] = 'SINGLE'


sra_need_redownload = []
check = True
for study_acc, p in tqdm(all_list.items()):
    l = f'{study_acc}/rawdata/srr_ids.list'
    if study_acc in removed_sra:
        continue
    if p == 'PAIRED':
        sra = download(l, pe=1, check=check)
        sra_need_redownload.append(sra)
    else:
        sra = download(l, check=check)
        sra_need_redownload.append(sra)
sra_need_redownload = [_ for _ in sra_need_redownload if _]

check = False
for study_acc in tqdm(sra_need_redownload):
    p = all_list[study_acc]
    l = f'{study_acc}/rawdata/srr_ids.list'
    if study_acc in removed_sra:
        continue
    if p == 'PAIRED':
        sra = download(l, pe=1, check=check)
        sra_need_redownload.append(sra)
    else:
        sra = download(l, check=check)
        sra_need_redownload.append(sra)

