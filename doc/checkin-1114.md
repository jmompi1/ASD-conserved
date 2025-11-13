# Check in #2 (11/24) 

**How you’ve addressed prior feedback (2 pts)**
Since our last submission we have decided to stick to analyzing only two raw single-cell RNA-seq datasets. We also are in the process of expanding our gene list to include 10 genes of both categories and we hope to identify trends in the expression of these genes across our three species. 

**New progress since last submission (4 pts)**
Since our last submission we have analyzed the raw Drosophila brain single-cell RNA seq data (Drosophila_single-cell_analysis.R), labelled the three primary cell types, and generated UMAPs of the expression of 5 of our original 6 genes. One gene was not in the annotation file the authors used when alligning their data so it doens't seem to be present in the data. 

Here is how we cleaned the data and determined the identities of the 19 clusters generated. 

**Quality control:** 
Selected cells where nFeature_RNA > 1500, nFeature_RNA < 4000 & the percent of mitochondrial genes was less than 10 %
## Pre-QC:
<img width="1222" height="816" alt="preQC" src="https://github.com/user-attachments/assets/3bac8170-8359-4be7-b75b-bdd0fad05c8e" />

## Post-QC:
<img width="1222" height="816" alt="postQC" src="https://github.com/user-attachments/assets/0942c575-4a93-4c6a-98ce-50fa811b0e5e" />

**Cluster labelling:** 
## Neural markers: 
<img width="1222" height="816" alt="elav" src="https://github.com/user-attachments/assets/9e8e7e84-cfcf-46c3-a3e7-736fe883ad39" />
<img width="1222" height="816" alt="ChAT" src="https://github.com/user-attachments/assets/66c0e9f3-ad7e-4156-8c69-de4303004084" />

## Glial markers: 
<img width="1222" height="816" alt="repo" src="https://github.com/user-attachments/assets/4d0f95e0-238c-4b4b-824e-2eaa87158417" />
<img width="1222" height="816" alt="loco" src="https://github.com/user-attachments/assets/97b8e5f2-e680-4939-9cc4-f08423852a7c" />
<img width="1222" height="816" alt="gcm" src="https://github.com/user-attachments/assets/456f04b5-4555-46b7-98a0-62966c2912b1" />


## Neural progenitor markers: 
<img width="1222" height="816" alt="N" src="https://github.com/user-attachments/assets/34d56d48-3e8f-4d1e-bb38-c78188ecf18e" />
<img width="1222" height="816" alt="dpn" src="https://github.com/user-attachments/assets/3efb37b7-90a1-435d-bcaf-9d689dba265c" />

## Preliminary cluster labelling: 
<img width="1222" height="816" alt="prelim_cluster" src="https://github.com/user-attachments/assets/89ecce9c-51bb-4c33-ba3a-fbfd5e0401cd" />

**Expression patterns of 5/6 of the original genes:** 

## Nrx-1 (causal) 
<img width="1222" height="816" alt="nrx-1" src="https://github.com/user-attachments/assets/a514149b-d00e-4092-8776-1159339c7fc0" />

Nrx-1 seems to be pretty uniformly expressed in neurons of the developping L1 larval brain. 

## rg (causal)
<img width="1222" height="816" alt="rg" src="https://github.com/user-attachments/assets/cf81521d-3c0e-4c0d-a114-39506b76ecec" />

rg seems to be pretty lowly expressed uniformly in neurons at this developmental stage. 

## kis (causal)
<img width="1222" height="816" alt="kis" src="https://github.com/user-attachments/assets/897774d0-9b44-4c23-8c2e-12c1125ca165" />

kis doesn't seem to be preferentially expressed in the developping L1 larval brain. 

## vnd (related)
<img width="1222" height="816" alt="vnd" src="https://github.com/user-attachments/assets/b54db38c-257a-4fee-a2ce-69bebfdba900" />

vnd seems to be pretty lowly expressed uniformly in neurons at this developmental stage. 

## SK(related)
<img width="1222" height="816" alt="SK" src="https://github.com/user-attachments/assets/de50abda-2fb1-4a64-928e-4dfef6cda0b6" />

SK seems to be pretty uniformly expressed in neurons of the developping L1 larval brain. 

**Project Organization (2 pts)**
Moving forward we want to expand our list of genes of interest to 10 for each gene group. I also want to work on identifying neural subtypes in the Drosophila Single-cell RNA-seq data. We are also in process of analyzing the raw C. elegans single-cell RNA-seq data. 

**Struggles you are encountering and questions you would like advice on (2 pts)**
We would love to get some advice moving forward on how to compare the expression patterns between organisms, and some statistical tests we could look into.
