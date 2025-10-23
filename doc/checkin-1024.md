<img width="579" height="552" alt="Screenshot 2025-10-22 at 11 07 08 PM" src="https://github.com/user-attachments/assets/d86d7a4f-b8ec-4117-8c82-ab769c30e641" /># Check in #1 (10/24) 

**How you’ve addressed prior feedback (2 pts)**
Moving forward we have decided to pull the majority of our genes from one primary GWAS study (Grove J, et al., 2019) as per Rajiv's suggestion. We have also decided to use SFARI Gene as a consensus for gene scoring, which was another great ressource recommended to us by Rajiv. We also decided to limit ourselves to the full analysis of 2 model organisms brain development single cell data for now as we start to get familiar with single-cell analysis. Since a comprehensive single-cell RNA-seq atlas with a user-friendly interface is available online for zebrafish (https://zebrahub.sf.czbiohub.org/transcriptomics) we will only analyze the raw Drosophila and C.elegans sc-RNA seq data. 

**New progress since last submission (4 pts)**
Since our last submission we have decided to focus on the expression of 6 genes in all three model organisms, three genes with a SFARI gene score of 3 and three genes with a SFARI score of 1. We also made sure to only select genes with confident gene orthoglues in all three model organisms. As mentioned previously we primarily relied on the GWAS study (Grove J, et al., 2019) but we also pulled genes from a paper that analyzed rare de novo and inherited coding variants in 42,607 autism cases (Zhou X, et al., 2022). Genes with a SFARI score of 3 have a single reported de novo mutation likely causing protein truncation, come from unreplicated association studies, or are identified through rare inherited mutations without rigorous statistical analysis. While genes with a SFARI score of 1 are genes who's mutations have thoroughly been discribed as contributing to ASD. By selecting genes from both these subsets we hope to describe some preliminary trends in expression patterns between causal genes (score of 1) and related genes (score of 3) in our three model organisms of interest during brain development. Here are the genes we selected:

**Causal genes:** 
-  **NRXN1**: nrx-1 (C. elegans), Nrx-1 (Drosophila), nrxn1a & nrx1b (Zebrafish)
-  **NBEA**: sel-2 (C. elegans), rg (Drosophila), nbeab (Zebrafish)
-  **CHD8***: chd-7 (C. elegans), kis (Drosophila), chd8 (Zebrafish)

**Related genes:**
- **NKX2-2**: ceh-22 (C. elegans), vnd (Drosophila), nkx2.2a (Zebrafish)
- **KCNN2**: kcnl-2 (C. Elegans), SK (Drosophila), kcnn2 (Zebrafish)
- **YWHAZ**: ftt-2 & par-5 (C. Elegans), 14-3-3ζ (Drosophila), ywhaz (Zebrafish)

In the meantime, we have also selected a Drosophila (Brunet Avalos C et al., 2019) and a C. elegans (Cao J et al., 2017) single-cell RNA-seq dataset, which we have began to process the raw data from. The Drosophila dataset is from first-instar larvae and the C. elegans dataset is from L2 larvae which are both stages that are crucial periods for brain development. However, in the meantime we have started to look at the expression profiles of our 6 genes in the zebrafish single-cell RNA-seq atlas a 5dpf, which marks the a key step in the development of neuronal activity. 

Preliminary analysis of the expression pattern of our genes in the Zebrafish sc-RNA seq expression atlas: 

# Cell type annotations of clusters where our genes of interest seem to be expressed: 
<img width="729" height="475" alt="Screenshot 2025-10-22 at 10 56 57 PM" src="https://github.com/user-attachments/assets/e2988539-8409-425d-80db-01e6bebcc27f" />

# nrx-1a (causal) 
<img width="514" height="474" alt="Screenshot 2025-10-22 at 11 00 50 PM" src="https://github.com/user-attachments/assets/acc63722-b136-4ba1-95d6-1b01f86dfdee" />
nrx-1 seems to primarily be expressed in the olfactory bulb, the optic tectum, in neurons (including retinal neurons), and the ventral forebrain. 

# nbeab (causal)
<img width="615" height="552" alt="Screenshot 2025-10-22 at 11 05 27 PM" src="https://github.com/user-attachments/assets/668e2e2a-0dd5-465f-867a-65dd99ac2ef1" />
nbeab seems to be lowly expressed at this developmental stage.

# chd8 (causal)
<img width="595" height="556" alt="Screenshot 2025-10-22 at 11 06 23 PM" src="https://github.com/user-attachments/assets/f32f01de-8c78-45e9-af84-fb3e25e38c9a" />
chd8 doesn't seem to be preferentially expressed in any cell type at this stage in development. 

# nkx2.2a (related)
<img width="586" height="550" alt="Screenshot 2025-10-22 at 11 07 26 PM" src="https://github.com/user-attachments/assets/dce35673-35e7-45fd-bc1e-0a908ddf0119" />
nkx2.2a seems to be preferentially expressed in the hindbrain, the diecephalon and in radial glia. 

# kcnn2 (related)
<img width="549" height="556" alt="Screenshot 2025-10-22 at 11 11 15 PM" src="https://github.com/user-attachments/assets/a5cdc685-d883-41e6-84fb-5f95efd309ac" />
kcn22 seems to be lowly expressed at this developmental stage.

# ywhaz (related)
<img width="556" height="554" alt="Screenshot 2025-10-22 at 11 11 39 PM" src="https://github.com/user-attachments/assets/cf413b71-dbd8-4bcf-ad2b-f96af417eb33" />
ywhaz seems to be highly expressed in the ventral forebrain, in neurons (including retinal), and in forebrain gabaergic neurons. 


**Project Organization (2 pts)**
*Struggles you are encountering and questions you would like advice on (2 pts)*
We would love to get some advice moving foreward on how we should go about analyzing our single-cell data and would love any aditional ressources or available datasets that we might be overlooking. 
