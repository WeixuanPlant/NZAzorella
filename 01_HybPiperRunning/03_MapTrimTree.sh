hybpiper retrieve_sequences --targetfile_dna mega353.fasta supercontig --sample_names namelist.txt --fasta_dir SupercontigsOutput
grep ">" -c *_supercontig.fasta | sed 's/_supercontig.fasta:/ /' | sort -g -k 2 | awk '$2 >= 15  { print $1}' > aln132_n15.txt

for i in $(cat aln132_n15.txt); do  awk '{ print $1}' "${i}"_supercontig.fasta | \
 awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | \
 sed '/^>/! s/N//g' > "${i}"_supercontig_fixed.fasta; done


#https://www.cyberciti.biz/faq/how-to-remove-carriage-return-in-linux-or-unix/
sed 's/\r$//' Treenamen132.csv > Treenamen132_fixed.csv
 
 
 
 
#!/bin/bash -e
#SBATCH --job-name=angiosperms353_aln  #job name (shows up in the queue)
#SBATCH --account=massey02696
#SBATCH --cpus-per-task=10   # number of CPUs per task (1 by default)
#SBATCH --time=10:00:00 #Walltime (HH:MM:SS)
#SBATCH --mem=10G  #Memory in GB 

 
#align the individual gene sequences with maff
module load MAFFT/7.429-gimkl-2020a
for word in $(cat aln132_n15.txt); do mafft --auto --thread 10 "${word}"_supercontig_fixed.fasta > "${word}"_aln.fasta; done


module load trimAl/1.4.1-GCC-9.2.0
for word in $(cat aln132_n15.txt); do echo "${word}"; trimal -in "${word}"_aln.fasta -out "${word}"_trim.fasta -gt 0.7; done 

grep -c '>' *_trim.fasta | sed 's/_trim.fasta:/ /' | sort -g -k 2 | awk '$2 >= 0  { print $1, $2}'

for word in $(cat ../aln132_n15.txt)
do
echo "${word}"
sed 's/,/\t/g' ../Treenamen132_fixed.csv | while read a b; do sed -i "s/$a/$b/" "${word}"_trim.fasta; done
done


for f in *.fastq.gz; do 
mv "$f" "$(echo "$f" | sed 's/fastq/fq/')"; done


while read i; do echo $i |grep -o "nomerge" | wc -l;  done < loci_gt07_n200_myoso.treefile

module load IQ-TREE/2.2.0.5-gimpi-2022a



##########################################################################
##########################################################################################################
##########################################################################################################
