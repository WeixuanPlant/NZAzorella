# Phylogenomics analysis of New Zealand *Azorella* 
Resolving Reticulate Evolutionary Histories of Polyploid Species of *Azorella* (Apiaceae) Endemic to New Zealand 

Manuscript link: to be updated

<p float="left">
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAcyno.jpg" height="150" />
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAcyno2.jpg" height="150" /> 
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAroughii.jpg" height="150" /> 
  <img src="https://github.com/WeixuanPlant/NZAzorella/blob/main/Supplymentary/filedAroughii2.jpg" height="150" /> 
</p>

#
### Abstract 
Genera with species of multiple ploidal levels provide models to understand successive rounds of whole genome duplication leading to intricate reticulate relationships of polyploid plant species. Here, we studied 17 polyploid taxa (species, subspecies, or varieties) in *Azorella* sections *Schizeilema* and *Stilbocarpa* (Apiaceae) that are mostly endemic to New Zealand. Species in these groups have diverse habits and leaf morphologies, distinct distributional ranges, and varying ploidal levels (4x, 6x, and 10x). Our goals were to determine the origins of the higher-level polyploids (6x and 10x), resolve species relationships, and assess the biogeography of the New Zealand *Azorella* species. We reconstructed a species tree from Hyb-Seq Angiosperms353 baits data for 131 individuals of all 17 taxa. We also reconstructed nrDNA and plastome phylogenies for 105 individuals using genome skim sequencing. Phylogenomic analysis of Hyb-Seq data using locus-based species trees and SNP-based approaches, together with comparison of nrDNA and plastome trees, showed that species diversification within New Zealand may relate to multiple origins from South America, which has been further shaped by additional rounds of polyploidy as well as hybridization or introgression. The two *Azorella* sections in New Zealand likely resulted from different biogeographic events from South America - one to the subantarctic islands (section Stilbocarpa) and a second to the South Island (section Schizeilema). In addition, within section *Schizeilema*, species have dispersed from the South Island (New Zealand) to Australia, the subantarctic islands, and the North Island (New Zealand). Our combined approach of phylogenomic analyses of plastome and nuclear locus-based data, together with SNP-based network approaches allowed us to determine the origins of some higher-level polyploids in New Zealand *Azorella* and revealed a more complex picture of historical and ongoing polyploidy and hybridization within these lineages. 

#
### Gene tree based  phylogeny inference 

We extracted Angiosperm353 loci via [HybPiper](https://github.com/mossmatters/HybPiper) using target-enrichment sequencing, and high-copy marker of plastome and nrDNA were extracted via [Getorgenalla](https://github.com/Kinggerm/GetOrganelle). All bioinformatic were performed via New Zealand eScience Infrastructure [NeSI](https://www.nesi.org.nz/)
 
### A353 loci sequence extraction and gene trees reconstruction
#### Raw reads trimming 
```
module load Trimmomatic/0.39-Java-1.8.0_144

for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/01_trimmed/*_1.fq.gz
	do
	withpath="${file}"
	filename=${withpath##*/}
	base="${filename%**_1.fq.gz}"
	echo "${base}"
	trimmomatic PE -threads 25 "${base}"_1.fq.gz "${base}"_2.fq.gz \
	"${base}"_1P.fq.gz "${base}"_1UP.fq.gz \
	"${base}"_2P.fq.gz "${base}"_2UP.fq.gz \
	ILLUMINACLIP:TruSeq3-PE-2.fa:2:20:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done
```

#### Supercontigs retreving via HybPiper
```
module load Python/3.9.9-gimkl-2020a
module load HybPiper/2.0.1rc-Miniconda3
module load Parallel/20200522

for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/02_pairedtrimed/*_1P.fq.gz
	do
	withpath="${file}"
	filename=${withpath##*/}
	base="${filename%*_1P.fq.gz}" 
	echo "${base}"
	
	hybpiper assemble --run_intronerate \
	--readfiles /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/02_pairedtrimed/"${base}"*P.fq.gz \
	--targetfile_dna mega353.fasta --bwa \
	--prefix "${base}"_nomerge  --no_padding_supercontigs \
	--timeout_assemble 600 --paralog_min_length_percentage 0.5
	done
```
### High copy gene sequence extraction
#### Raw reads trimming
```
module load Trimmomatic/0.39-Java-1.8.0_144

for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_GenomeSkimming_NCBI_upload/01_trimmedreads/*_1.fq.gz
        do
        withpath="${file}"
        filename=${withpath##*/}
        base="${filename%**_1.fq.gz}"
        echo "${base}"
        trimmomatic PE -threads 25 "${base}"_1.fq.gz "${base}"_2.fq.gz \
        "${base}"_1P.fq.gz "${base}"_1UP.fq.gz \
        "${base}"_2P.fq.gz "${base}"_2UP.fq.gz \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:20:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done
```
#### Retreiving high-copy markers via Getorganelle 
```
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /scale_wlg_persistent/filesets/project/massey02696/biopython
cd /nesi/nobackup/massey02696/WeixuanData/Azorella_GenomeSkimming_NCBI_upload/03_getorgenlle


for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_GenomeSkimming_NCBI_upload/02_trimmedreads/*_1P.fq.gz 
	do 
	withpath="${file}" 
	filename=${withpath##*/} 
	base="${filename%*_1P.fq.gz}" 
	echo "${base}" 
	
	get_organelle_from_reads.py -1 ../02_trimmedreads/"${base}"_1P.fq.gz -2 ../02_trimmedreads/"${base}"_2P.fq.gz -o plastome_output/"${base}" -R 15 -k 21,45,65,85,105 -F embplant_pt
	get_organelle_from_reads.py -1 ../02_trimmedreads/"${base}"_1P.fq.gz -2 ../02_trimmedreads/"${base}"_2P.fq.gz -o nrdna_output/"${base}" -R 10 -k 35,85,115 -F embplant_nr
	done
```
#
### SNPs based phylogeny analysis 
#### SNP calling via remapping the reads back to the same reference for 23 representative individuals
```
#The reference fasta gene names need to be changed according to file names
#The script variantcall need to change fastaq to fq
#The GenotypetoPCA need to change the script expression & to || for GATK
#The plink needs to update the version 

module rest
module load GATK/4.1.4.1-gimkl-2018b
module load PLINK/1.09b6.16
module load SAMtools/1.8-gimkl-2018b
module load BWA/0.7.17-gimkl-2017a
module load BCFtools/1.9-GCC-7.4.0
module load Python/3.7.3-gimkl-2018b
module load WhatsHap/1.1-gimkl-2020a

bash variantcall.sh Azho-AK16_S7_L001.supercontigs.fasta Azho-AK16_S7_L001 #please note that sample name do not have the slash / at the end
bash GenotypesToPCA.sh Azho-AK16_S7_L001.supercontigs.fasta Azho-AK16_S7_L001
bash plink_stats.sh Azho-AK16_S7_L001
```

#
### This pipeline still underconstruction!

