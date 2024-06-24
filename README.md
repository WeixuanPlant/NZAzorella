# Phylogenomics analysis of New Zealand *Azorella* 
Hyb-Seq (Angiosperms353 baitset) and Genome-skimming sequencing of New Zealand *Azorella*   

Manuscript see: LINK inserted here (https://www.biorxiv.org/)

<p float="left">
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAcyno.jpg" height="150" />
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAcyno2.jpg" height="150" /> 
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAroughii.jpg" height="150" /> 
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAroughii2.jpg" height="150" /> 
</p>

### Gene tree based  phylogeny inference 

#### Single copy gene sequence extraction
▪️ Pop gene analysis (PCA, NJ-tree, Structure)
1. All AD1 65 samples' raw reads were trimmed with [Trimmomatic](https://github.com/usadellab/Trimmomatic.git).
3. Following [sentieon-dnaseq](https://github.com/Sentieon/sentieon-dnaseq.git), trimmmed reads mapping; GVCF calling; VCF calling; SNPs filtering.
4. Using biallelic SNPs to estimate population genetic groups via [PLINK](https://www.cog-genomics.org/plink/) (PCA) and [LEA](https://bioconductor.org/packages/release/bioc/html/LEA.html).
5. Including additional two AD4 samples as outgroup, and calling bialleic SNPs from the 'combined' VCF with 65 AD1 samples, which include variable and invariable sites.

#### High copy gene sequence extraction


### SNPs based phylogeny analysis 

#### SNPs extraction from Hyb-Seq data

#### SNPs extraction from Hyb-Seq data



### Please cite the paper: 

