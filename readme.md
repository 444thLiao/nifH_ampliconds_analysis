# Brief instructions of nifH amplicons analysis


All path implemented in scripts should be carefully modified accordingly.


### 1. download SRR (SRA runs) 

use the script `SRR_download.py` to download retrieved metadata data.


### 2. run amplicons profilling pipelines

For amplicons data with quality information, we applied qiim2 pipleines to process it . (which implemented in `workflows/ampliconds_PE_analysis.py` and `workflows/ampliconds_data_process.py`.


But for amplicons data without quality information, we applied `vsearch` to profilling those sequencing data.


see script `run_ampliconds_profilling.py`

### 3. phylogenetic placements

Leveraging to the generated representative sequences from profilling pipelines, we applied the alignment-based phylogenetic placements (**EPA-ng**)


see script `phylogenetic_placements.py`


### 4. visualizations

After phylogenetic placements, the generated `jplace` files are used to classification. You coul see the several functions (such as `classification_criteria`) at script`visualizations.py`


Furthermore, you could also refer to the ipynb named `nifH visualizations.ipyn`



Reference and detailed methods see publications 

[Tao J, Wang S, Liao T, et al. Evolutionary origin and ecological implication of a unique nif island in free-living Bradyrhizobium lineages[J]. The ISME Journal, 2021: 1-12.](https://doi.org/10.1038/s41396-021-01002-z)


For More information

Please contact: l0404th@gamil.com
