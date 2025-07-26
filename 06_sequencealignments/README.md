## Final aligned sequences 

### Target-enriched trimmed gene alignments for selected 345 individuals: [aln_n345.7z](https://github.com/WeixuanPlant/NZAzorella/blob/main/06_sequencealignments/aln_loci219_n23.7z) 

### Target-enriched trimmed gene alignments for selected 23 individuals: [aln_loci219_n23.7z](https://github.com/WeixuanPlant/NZAzorella/blob/main/06_sequencealignments/aln_n345.7z)

### raw alignment for nrDNA sequences: [nrDNA_sorted_aln_renamed.fasta](https://github.com/WeixuanPlant/NZAzorella/blob/main/06_sequencealignments/nrDNA_sorted_aln_renamed.fasta)

### raw alignment for plastome sequences: [plastome_n99complete_n6map_aln.fasta](https://github.com/WeixuanPlant/NZAzorella/blob/main/06_sequencealignments/plastome_n99complete_n6map_aln.fasta)


## Final reconstructed phylogenies

### Final nrDNA tree for 105 *Azorella* samples: [Azorella_nrDNA_sorted_aln_trimgt07_rerooted.treefile](https://github.com/WeixuanPlant/NZAzorella/blob/main/06_sequencealignments/Azorella_nrDNA_sorted_aln_trimgt07_rerooted.treefile)

### Final plastome tree for 105 *Azorella* samples: [Azorella_plastome_n99complete_n6map_aln_trim2_rerooted.treefile](https://github.com/WeixuanPlant/NZAzorella/blob/main/06_sequencealignments/Azorella_plastome_n99complete_n6map_aln_trim2_rerooted.treefile)

### Target-enriched gene trees: [loci_gt07_n15_gene345_bs10_reroot.treefile](https://github.com/WeixuanPlant/NZAzorella/blob/main/01_HybPiperRunning/loci_gt07_n15_gene345_bs10_reroot.treefile)

### Target-enriched species ASTRAL trees: [loci_gt07_n15_gene345_bs10_reroot.treefile](https://github.com/WeixuanPlant/NZAzorella/blob/main/01_HybPiperRunning/loci_gt07_n15_gene345_bs10_reroot_astral_reroot.treefile)

```
sed 's/,/\t/g' rename.csv | while read a b; do sed -i "s/$a/$b/" plastome_n99complete_n6map_aln.fasta ; done

seqkit mutate -w 0 plastome_n99complete_n6map_aln.fasta | sed 's/-//g'  > plastome_n99complete_n6map_ncbi.fasta

awk '/^>/ {
    orig_header = $0;
    split(substr($0, 2), parts, "_");
    prefix = parts[1];
    if (prefix ~ /^Az/) {
        species = tolower(substr(prefix, 3));
        organism = "Azorella " species;
    } else {
        organism = prefix;
    }
    print orig_header "  [organism=" organism "]  [mol_type=Chloroplast DNA]";
    next
}
{ print }'  plastome_n99complete_n6map_ncbi.fasta > plastome_n99complete_n6map_ncbi2.fasta

sed -i 's/Azorella haastii/Azorella haastii subsp. haastii/' plastome_n99complete_n6map_ncbi2.fasta 
sed -i 's/Azorella cyanopetala/Azorella haastii subsp. cyanopetala/' plastome_n99complete_n6map_ncbi2.fasta

table2asn -i plastome_n99complete_n6map_ncbi2.fasta


```
