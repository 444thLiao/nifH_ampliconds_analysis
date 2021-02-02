"""
Perform the amplicons profilling analysis
"""
import multiprocessing as mp
from os.path import *
from tqdm import tqdm
import os

import warnings
warnings.filterwarnings("ignore")

from .SRR_download import check_done,all_list,removed_sra,no_quality_ids
from .workflows.amplicons_data_process import se_workflow
from .workflows.amplicons_PE_analysis import pe_workflow

def check_wrong_layout(infile, pe=False):
    srrids = open(infile).read().split('\n')
    srrids = [_ for _ in srrids if _]
    sra_id = infile.split('/')[-3]

    wrong_layout = False
    for srr in srrids:
        if pe:
            if exists(f"{dirname(infile)}/{srr}_1.fastq") and (not exists(f"{dirname(infile)}/{srr}_2.fastq")):
                os.system(
                    f"mv {dirname(infile)}/{srr}_1.fastq {dirname(infile)}/{srr}.fastq")
                wrong_layout = True
            elif exists(f"{dirname(infile)}/{srr}.fastq"):
                wrong_layout = True
    if wrong_layout:
        tqdm.write(f"{sra_id} should be Single")
        return False
    return True


wrong_sra_ids = []
for sra_id, p in all_list.items():
    l = f'/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/{sra_id}/rawdata/srr_ids.list'
    if sra_id in removed_sra:
        continue
    if p == 'PAIRED':
        if not check_wrong_layout(l, pe=1):
            wrong_sra_ids.append(sra_id)


def run(args):
    unit_run(*args)


def unit_run(sra_id, pe=False, quick=False):
    if pe:
        pe_workflow(sra_id, quick=quick)
    else:
        se_workflow(sra_id, quick=quick)


dry_run = False
processed_sra_ids = []
for sra_id, p in all_list.items():
    l = f'/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/{sra_id}/rawdata/srr_ids.list'
    all_srr = open(l).read().strip('\n').split('\n')
    if sra_id in removed_sra:
        continue
    if p == 'PAIRED':
        if check_done(l, pe=1):
            tqdm.write('\t'.join([sra_id, p, str(len(all_srr))]))
            if not dry_run:
                processed_sra_ids.append((sra_id, True))
                # os.chdir(dirname(dirname(l)))
                # pe_workflow(sra_id)
    else:
        if check_done(l):
            tqdm.write('\t'.join([sra_id, p, str(len(all_srr))]))
            if not dry_run:
                processed_sra_ids.append((sra_id, False))
                # os.chdir(dirname(dirname(l)))
                # se_workflow(sra_id)

with mp.Pool(processes=3) as tp:
    r = list(tqdm(tp.imap(run, processed_sra_ids), total=len(processed_sra_ids)))


for sra in no_quality_ids:
    if all_list[sra] == 'SINGLE':
        se_workflow(sra, quality=False)
    else:
        pe_workflow(sra, quality=False)