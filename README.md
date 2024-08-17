# SISRS Amaranthus Whole Genome Barcode Discovery Pipeline

![Screenshot 2024-08-07 at 9 31 50 PM](https://github.com/user-attachments/assets/ceb749bf-acbe-4223-bcbe-f7808cb1c60b)

Referenced from the original SISRS : https://github.com/BobLiterman/2024_SISRS_General.git
Here's a structured README for your SISRS project, formatted similarly to the example you provided:

---

# Running SISRS for Amaranthus Genomic Barcode Development

### Lubaina Kothari   
University of Guelph - Masters of Bioinformatics Internship

This repository provides a comprehensive guide for using the SISRS pipeline to identify diagnostic SNPs for species-specific orthologs in the *Amaranthus* genus. It outlines each step from data organization to final SNP analysis and includes detailed commands and scripts.

---

## Dependencies

Before you begin, ensure you have the following software installed:

- **FastQC** (v0.11.9)
- **BBTools** (BBDuk v38.84)
- **GetOrganelle** (v1.7.5)
- **Ray** (v2.3.1)
- **Bowtie2** (v2.4.2)
- **Samtools** (v1.10+)
- **Python** (v3.7+)
- **IQ-TREE** (v2.0.6)
- **BEDTools** (v2.29.2)

---

## Step 1: Quality Filtering and Trimming

### Description
Perform quality control on raw sequencing data using FastQC, followed by trimming of adapters and low-quality bases using BBDuk.

### Command
```bash
# FastQC quality check
fastqc *.fastq

# Adapter trimming and quality filtering with BBDuk
bbduk.sh in1=R1.fastq in2=R2.fastq out1=R1_trimmed.fastq out2=R2_trimmed.fastq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=10 minlength=50
```

### Output
- Quality reports (`.html` and `.zip` files)
- Trimmed FASTQ files (`R1_trimmed.fastq`, `R2_trimmed.fastq`)

---

## Step 2: Pan Plastid Genome Assembly and Nuclear Read Extraction

### Description
Construct a pan-plastid genome to enrich nuclear reads by removing plastid contamination. Use GetOrganelle for assembly and BBDuk for extraction.

### Command
```bash
# Assemble plastid genomes using GetOrganelle
get_organelle_from_reads.py -1 R1_trimmed.fastq -2 R2_trimmed.fastq -o output_directory -R 10 -k 21,45,65,85,105

# Combine assembled plastids to form the pan-plastid genome
cat output_directory/*.fasta > pan_plastid.fasta

# Extract plastid reads using BBDuk
bbduk.sh in1=R1_trimmed.fastq in2=R2_trimmed.fastq out1=R1_nuclear.fastq out2=R2_nuclear.fastq ref=pan_plastid.fasta k=31 hdist=1
```

### Output
- Pan-plastid genome (`pan_plastid.fasta`)
- Nuclear-enriched FASTQ files (`R1_nuclear.fastq`, `R2_nuclear.fastq`)

---

## Step 3: Composite Genome Assembly

### Description
Organize reads by species and subset them equally to ensure balanced representation in the composite genome. Use the Ray assembler to construct the composite genome.

### Command
```bash
# Organize reads into species-specific folders
mkdir -p Species_A; mv R1_nuclear.fastq R2_nuclear.fastq Species_A/
mkdir -p Species_B; mv R1_nuclear.fastq R2_nuclear.fastq Species_B/

# Subset reads using SISRS read_subsetter.py script
python read_subsetter.py -d Species_A -s 10 -g 500M
python read_subsetter.py -d Species_B -s 10 -g 500M

# Assemble the composite genome with Ray
mpirun -np 4 Ray -k 31 -o composite_genome -s Species_A/subset_R1.fastq -p Species_A/subset_R1.fastq Species_A/subset_R2.fastq
```

### Output
- Composite genome assembly (`composite_genome/Contigs.fasta`)

---

## Step 4: Mapping and Indexing

### Description
Map the sample reads to the composite genome using Bowtie2 and create an index to facilitate efficient mapping.

### Command
```bash
# Build Bowtie2 index
bowtie2-build composite_genome/Contigs.fasta composite_genome

# Map reads back to the composite genome
bowtie2 -x composite_genome -1 R1_nuclear.fastq -2 R2_nuclear.fastq -S sample_mapped.sam

# Convert SAM to BAM and sort
samtools view -bS sample_mapped.sam | samtools sort -o sample_mapped_sorted.bam
```

### Output
- Indexed composite genome (`.bt2` files)
- Sorted BAM files (`sample_mapped_sorted.bam`)

---

## Step 5: Sample-Specific Ortholog Generation

### Description
Generate orthologous sequences for each species by mapping, identifying fixed alleles, and filtering based on missing data and gaps.

### Command
```bash
# Generate fixed allele alignment
python Output_SISRS.py -a composite_genome/Contigs.fasta -b Species_A.bam -c Species_B.bam -o fixed_alleles_alignment

# Filter alignment for missing data and gaps
python filter_nexus_for_missing.py fixed_alleles_alignment.nex 0
```

### Output
- Filtered alignments (`filtered_alignment.nex`, `filtered_alignment.phylip-relaxed`)

---

## Step 6: Phylogenetic Tree Construction

### Description
Construct phylogenetic trees using the filtered alignments with IQ-TREE.

### Command
```bash
iqtree -nt AUTO -s filtered_alignment.phylip-relaxed -m MFP+MERGE+ASC -bb 1000
```

### Output
- Phylogenetic tree files (`.treefile`, `.contree`, `.log`)

---

## Step 7: Species-Specific SNP Identification

### Description
Identify species-specific SNPs using the filtered alignments and generate a BED file for downstream analysis.

### Command
```bash
# Extract SNPs using alignment_slicer.py
python alignment_slicer.py filtered_alignment.phylip-relaxed species_specific_SNPs.bed

# Convert SNPs to BED format
awk '{print $1 "\t" $2-1 "\t" $2}' species_specific_SNPs.bed > species_specific_SNPs.bed
```

### Output
- BED file with species-specific SNPs (`species_specific_SNPs.bed`)

---

## Step 8: Downstream Analysis and Validation

### Description
Use BEDTools to analyze dense SNP regions and Rboretum for SNP validation and species classification.

### Command
```bash
# Identify dense SNP regions with BEDTools
bedtools cluster -i species_specific_SNPs.bed

# Perform SNP validation using Rboretum in R
Rscript validate_SNPs.R species_specific_SNPs.bed
```

### Output
- SNP clusters (`.bed`)
- Validation results (`SNP_validation_report.txt`)

---

This README provides a detailed guide to execute each step of the SISRS pipeline, including the necessary commands and expected outputs. For more information and additional scripts, please refer to the scripts in this repository.# 1. Quality Control of Sequences

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
