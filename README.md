# SISRS Amaranthus Whole Genome Barcode Discovery Pipeline

![Screenshot 2024-08-07 at 9 31 50 PM](https://github.com/user-attachments/assets/ceb749bf-acbe-4223-bcbe-f7808cb1c60b)

# 1. Quality Control of Sequences

1.1 Perform FastQC Check Run FastQC to check the quality of the raw sequencing data. 
mkdir -p FastQC_Reports fastqc -o FastQC_Reports raw_data/*.fastq.gz

1.2 Adapter Trimming with BBDuk Trim adapters and perform light quality trimming using BBDuk. 
bbduk.sh in1=raw_data/sample_R1.fastq.gz in2=raw_data/sample_R2.fastq.gz \ out1=trimmed/sample_trimmed_R1.fastq.gz out2=trimmed/sample_trimmed_R2.fastq.gz \ ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

# 2. Nuclear Enrichment (Removal of plastid reads)

2.1 Run GetOrganelle Assemble the pan plastid genome using GetOrganelle. 
get_organelle_from_reads.py -1 trimmed/sample_trimmed_R1.fastq.gz -2 trimmed/sample_trimmed_R2.fastq.gz \ -o getorganelle_output -F embplant_pt -R 10 --cpu 8

Concatenate files for a final pan-plastid assembly.

2.2 Extract Pan Plastid Reads Use BBDuk to Extract Pan Plastid Reads Extract reads matching the pan plastid genome from each sample. 
bbduk.sh in1=trimmed/sample_trimmed_R1.fastq.gz in2=trimmed/sample_trimmed_R2.fastq.gz \ out1=pan_plastid/sample_plastid_R1.fastq.gz out2=pan_plastid/sample_plastid_R2.fastq.gz \ ref=pan_plastid_genome.fa k=31 hdist=1

# 3. SISRS Pipeline

3.1 Organize Samples into Taxonomic Folders Organize your samples into suspected taxonomic folders. 
mkdir -p SISRS_Project/{Species1,Species2,...} mv sample1.fastq.gz SISRS_Project/Species1/ mv sample2.fastq.gz SISRS_Project/Species2/

3.2 Subsetting and Composite Genome Assembly Subset reads equally among taxa and assemble the composite genome. Specify genome size and target coverage. 10x is advisable and works well for skimmed data. python Read_Subsetter.py -g 4500000000 -c 10 -r Input -o Output Note: Read_Subsetter script expects '.fastq.gz' with the naming convention '_1.fastq.gz/_2.fastq.gz' for forward and reverse reads respectively. If this is not the case for your files, please adjust the script using the below example: # If your input reads are _R1.fastq/_R2.fastq, you could run: sed -i 's/_1.fastq.gz/_R1.fastq/g' Read_Subsetter.py sed -i 's/_2.fastq.gz/_R2.fastq/g' Read_Subsetter.py sed -i 's/fastq.gz/fastq/g' Read_Subsetter.py

3.3 Assemble Composite Genome Using Ray Use Ray for genome assembly. 
Ray -p SISRS_Project/Composite_reads -o Composite_genome Note: Ray output folders must start with 'Ray_'

3.4 Prepare Composite Genome Use script 'Prepare_Composite.sh' which calls Genome_SiteLengths.py Contains the following order of processing: 
rename.sh in=SISRS_Project/Composite_Genome/Ray_Dir/Contigs.fasta out=SISRS_Project/Composite_Genome/Ray_Dir/Composite_Genome/contigs.fa prefix=SISRS addprefix=t trd=t

bowtie2-build contigs.fa contigs -p PROCESSOR_COUNT bbmap.sh ref=contigs.fa samtools faidx contigs.fa

python scripts/Genome_SiteLengths.py SISRS_Project/Composite_genome
