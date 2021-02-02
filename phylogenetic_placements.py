
from os.path import *
import os
from glob import glob
from tqdm import tqdm


ref_newick = realpath('./tree/iqtree.renamed.newick')
ref_phy = realpath('./ref.phy')
fa_iter_pattern = '/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/*/*_repr_dada2.fa'

cmds = []
for infa in tqdm(glob(fa_iter_pattern)):
    srr_id = basename(dirname(infa))
    infa = realpath(infa)
    
    odir_aln = f'./aln_f/{srr_id}' 
    if not exists(odir_aln):
        os.system(f"mkdir -p {odir_aln}")
    cmd1 = f"cd {odir_aln}; papara -t {ref_newick} -s {ref_phy} -q {infa} -r -n aln -j 30"
    if not exists(f"{odir_aln}/papara_alignment.aln"):
        cmds.append(cmd1)
    
    cmd2 = f"cd {odir_aln}; epa-ng --split {ref_phy} papara_alignment.aln"
    if exists(f"{odir_aln}/papara_alignment.aln") and (not exists(f"{odir_aln}/reference.fasta") ):
        cmds.append(cmd2)
    
    epa_output_dir = f"./epa_o/{srr_id}"
    if not exists(epa_output_dir):
        os.system(f"mkdir -p {epa_output_dir}")
    cmd3 = f"epa-ng --ref-msa {odir_aln}/reference.fasta --tree {ref_newick} --query {odir_aln}/query.fasta -T 30 -w {epa_output_dir} --model " + "GTR{0.7/1.8/1.2/0.6/3.0/1.0}+F+R10 --redo"
    if exists(f"{odir_aln}/reference.fasta") and (not exists(f"{epa_output_dir}/epa_result.jplace") ):
        cmds.append(cmd3)
        
        
