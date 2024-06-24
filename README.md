# Phylogenomics analysis of New Zealand *Azorella* 
Hyb-Seq (Angiosperms353 baitset) and Genome-skimming sequencing of New Zealand *Azorella*   

Manuscript see: LINK inserted here (https://www.biorxiv.org/)

<p float="left">
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAcyno.jpg" height="150" />
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAcyno2.jpg" height="150" /> 
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAroughii.jpg" height="150" /> 
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAroughii2.jpg" height="150" /> 
</p>

### Gene-based tree analysis 
▪️ Pop gene analysis (PCA, NJ-tree, Structure)
1. All AD1 65 samples' raw reads were trimmed with [Trimmomatic](https://github.com/usadellab/Trimmomatic.git).
   > trimmomatic PE -threads $thr $file1 $file2 $tDir/$name.R1.fq.gz $tDir/$name.U1.fq.gz $tDir/$name.R2.fq.gz $tDir/$name.U2.fq.gz ILLUMINACLIP:Adapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:75

3. Following [sentieon-dnaseq](https://github.com/Sentieon/sentieon-dnaseq.git), trimmmed reads mapping; GVCF calling; VCF calling; SNPs filtering.
   > vcftools --vcf $output2.vcf --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 1 --maf 0.05 --recode --recode-INFO-all --out $output2.variant

4. Using biallelic SNPs to estimate population genetic groups via [PLINK](https://www.cog-genomics.org/plink/) (PCA) and [LEA](https://bioconductor.org/packages/release/bioc/html/LEA.html).
5. Including additional two AD4 samples as outgroup, and calling bialleic SNPs from the 'combined' VCF with 65 AD1 samples, which include variable and invariable sites.

####  ▪️ Genetic variation comparison (Pi, Dxy, Fst, He, Fis, LD)
1. [Pixy](https://github.com/ksamuk/pixy.git) was applied to 


####  ▪️ Novel SNPs tabulating
1. Bcftools

####  ▪️ MK cotton population demographic analysis (PCA, Tajima's D, SFS, Ne)


### Please cite the paper: 

