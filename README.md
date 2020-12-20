# Andreev_Zhurbenko_fall_2020
 Assembly and analysis of organelle genomes of Iris species and Picea abies

Project goal: to write a pipeline for organelle genome assembly based on the reference. Assemble and analyze the mitochondrial genome of Picea abies and chloroplast genomes of irises (own sequence data)

Project tasks:
1.	Genome assembly of irises and Picea using reference;
2.	Phylogeny analysis of Iris species; 
3.	Search for nuclear mitochondrial DNA (NUMT) in Picea genome.


List of used tools:
1.	BWA (Version: 0.7.17-r1188);
2.	Blast (Version: 2.6.0+)
3.	datamash (Version: 1.7);
4.	IQ-TREE multicore (Version:  1.6.12)
5.	FastQC (Version: 0.11.9; Requirements: A suitable Java Runtime Environments);
6.	SAMtools. (Version: 1.11 (using htslib 1.11));
7.	SPAdes genome assembler (Version: 3.14.1 Requirements: 64-bit Linux system, 16 GB RAM, 50 GB free disc space).
8.	Quast (Version: 5.1.0rc1, e010ca46, Requirements: Linux (64-bit and 32-bit with slightly limited functionality) and macOS are supported; For the main pipeline: Python2 (2.5 or higher) or Python3 (3.3 or higher),  Perl 5.6.0 or higher, GCC 4.7 or higher, GNU make and ar, zlib development files; For the optional submodules: Time::HiRes perl module for GeneMark-ES (needed when using --gene-finding --eukaryote), Java 1.8 or later for GRIDSS (needed for SV detection) , R for GRIDSS (needed for SV detection))
9.  trimmomatic (Version: 0.39)

Iris genomes assembly

Work on server. Download the folder:
>scp -P 3622 -r -oHostKeyAlgorithms=+ssh-dss andreev@84.52.88.183:/mnt/mnt/andreevalexandr/Prject_Picea_Iris/irirs/trees/Iris_mafft.fas.treefile /home/pzhurb/Genomes/

Quality filtering using trimmomatic:
>java -jar trimmomatic-0.39.jar SE -phred33 ~/Iris_plastome/fastq/Fili_1_CTTGTA_L007_R1_001.fastq Fili_R1.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20

Running SPAdes on server:
>python /home/tools/SPAdes-3.14.1-Linux/bin/spades.py -s Fili_R1.fastq -o fili_fastq --isolate -t 14

Explore the read depth:
>samtools depth Fili_align_sort.bam | datamash sum 3 mean 3
mean – 110, max – 302

Index the reference with BWA:
>bwa index reference_coli.fna

Aligning on reference:
>bwa mem gat_ref.fasta contigs.fasta > alignments.sam; samtools view -S -b alignments.sam > alignments.bam; samtools sort alignments.bam -o alignments_sorted.bam; samtools index alignments_sorted.bam

Variant calling and creation of consensus alignment:
>bcftools mpileup -d 301 -L 301  -Ou -f gat_ref.fasta Fili_align_sort.bam | bcftools call -mv -Ou | bcftools filter -Oz -o calls1_ind_ins.vcf.gz | bcftools index calls1_ind_ins.vcf.gz

IQtree
>iqtree -s Iris_mafft.fas -bb 1000 -nt 18

Mitochondrial genome analysis
Picea abies nuclear genome 
> wget ftp://plantgenie.org/Data/ConGenIE//Picea_abies/v1.0/FASTA/GenomeAssemblies/Pabies1.0-genome.fa.gz

Picea abies mitochondrial genome 
>wget ftp://plantgenie.org/Publications/Sullivan2019//Picea_abies_mitochondrial_genome_annotation_curated.zip
>wget ftp://plantgenie.org/Publications/Sullivan2019//Picea_abies_mtDNA_assembly.zip
Unzip:
>unzip *.zip

Blast:
>makeblastdb -in ~/ Picea_abies_mtDNA_assembly.fa -dbtype nucl -parse_seqids (индексация);
>blastn -db Picea_abies_mtDNA_assembly.fa -quary Pabies1.0-genome.fa -out results.csv -outfmt "10 pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"

Searching for the contigs:
>quast.py -s -e NAD1.fas -r Pabies1.0-genome.fa -o Quast_Picea
>cat all_alignments_output.tsv |awk '{if($5!=unaligned)print $6 }' > id.txt
>grep -w -A 2 -f id.txt Pabies1.0-genome.fa --no-group-separator

