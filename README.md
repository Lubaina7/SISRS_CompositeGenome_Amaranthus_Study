# SISRS Amaranthus Whole Genome Barcode Discovery Pipeline

![Screenshot 2024-08-07 at 9 31 50 PM](https://github.com/user-attachments/assets/ceb749bf-acbe-4223-bcbe-f7808cb1c60b)

Referenced from the original SISRS : https://github.com/BobLiterman/2024_SISRS_General.git


---

# **Running SISRS for Amaranthus Genomic Barcode Development**

### Lubaina Kothari   
University of Guelph - Masters of Bioinformatics Internship  
In collaboration with the MIRL lab at CFIA

This repository provides a comprehensive guide for using the SISRS pipeline to identify diagnostic SNPs for species-specific orthologs in the *Amaranthus* genus. It outlines each step from data organization to final SNP analysis and includes detailed commands and scripts.

---

## **Set Up Environment and Dependencies**

Before you begin, ensure to have the following software installed. Using **mamba** is recommended to efficiently manage dependencies.

### Install Mamba
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh # By default, this installs in ~/miniforge3
```

### Create Your Environment
```bash
mamba create --prefix /path/to/env/location/ python=3.10
```

### Activate the Environment
```bash
mamba activate /path/to/env/location/
```

### Install Required Packages
```bash
mamba install -c conda-forge bioconda samtools bowtie2 bbmap fastqc trimmomatic
```

### Ensure the Following Packages are Installed:
- **FastQC** (v0.11.9)
- **BBTools** (BBDuk v38.84)
- **GetOrganelle** (v1.7.5)
- **Ray** (v2.3.1)
- **Bowtie2** (v2.4.2)
- **Samtools** (v1.10+)
- **Python** (v3.7+)
- **IQ-TREE** (v2.0.6)
- **BEDTools** (v2.29.2)

### Clone the SISRS Repository
Next, clone the SISRS repository to access all necessary scripts for the pipeline.
```bash
# Clone the SISRS repository
git clone https://github.com/BobLiterman/SISRS.git
cd SISRS
```

---

## **Step 1: Quality Filtering and Trimming**

### Description
Perform quality control on raw sequencing data using FastQC, followed by trimming of adapters and low-quality bases using BBDuk.

```bash
# FastQC quality check
fastqc *.fastq

# Adapter trimming and quality filtering with BBDuk
bbduk.sh in1=R1.fastq in2=R2.fastq out1=R1_trimmed.fastq out2=R2_trimmed.fastq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=10 minlength=50
```

### Output
- **Quality reports** (`.html` and `.zip` files)
- **Trimmed FASTQ files** (`R1_trimmed.fastq`, `R2_trimmed.fastq`)

---

## **Step 2: Running BBMerge for Plastid Assembly (Optional)**

### Description
Optionally, use BBMerge to merge paired-end reads, increasing the length of reads for better plastid assembly in GetOrganelle. Note: The reads can be left as paired and not merged as well for pan-plastid assembly.

```bash
# Loop through all R1 files to merge paired-end reads
for r1 in trimmed_*_R1_*.fastq.gz; do
  r2=${r1/_R1_/_R2_}
  merged_out="merged_${r1/_R1_001.fastq.gz/_merged.fastq.gz}"
  unmerged_out1="unmerged_${r1}"
  unmerged_out2="unmerged_${r2}"
  bbmerge.sh in1=$r1 in2=$r2 out=$merged_out outu1=$unmerged_out1 outu2=$unmerged_out2
done
```

---

## **Step 3: Pan Plastid Genome Assembly and Nuclear Read Extraction**

### Description
Construct a pan-plastid genome to enrich nuclear reads by removing plastid contamination. Use GetOrganelle for assembly and BBDuk for extraction.

```bash
# Activate the correct environment for GetOrganelle
source activate amaranthus_env

# Run GetOrganelle for each sample
OUTPUT_DIR="getOrganelle_output"
ORGANELLE_TYPE="embplant_pt"

mkdir -p $OUTPUT_DIR

for r1 in trimmed_*_R1_*.fastq.gz; do
  r2=${r1/_R1_/_R2_}
  get_organelle_from_reads.py -1 $r1 -2 $r2 -o $OUTPUT_DIR -k 21,45,65,85,105 -R 50 -F $ORGANELLE_TYPE --overwrite
done

# Combine assembled plastids to form the pan-plastid genome
cat $OUTPUT_DIR/*.fasta > pan_plastid.fasta

# Extract plastid reads using BBDuk
bbduk.sh in1=trimmed_R1.fastq.gz in2=trimmed_R2.fastq.gz outu1=R1_nuclear.fastq.gz outu2=R2_nuclear.fastq.gz ref=pan_plastid.fasta k=31 hdist=1 outu=R1_nuclear.fastq.gz
```

### Output
- **Pan-plastid genome** (`pan_plastid.fasta`)
- **Nuclear-enriched FASTQ files** (`R1_nuclear.fastq.gz`, `R2_nuclear.fastq.gz`)

---

## **Step 4: Subset Reads**

### Description
Organize reads by species and subset them equally to ensure balanced representation in the composite genome. Subsetting is performed based on desired coverage for the composite genome and the average genome size of the species. In this study, 10X coverage was used for a 0.5Gb average genome size in *Amaranthus* species. The custom SISRS script `read_subsetter.py` is used for this step.

### Note
If your input reads are not named `_1/_2.fastq.gz`, adjust `Read_Subsetter.py`:
```bash
sed -i 's/_1.fastq.gz/_R1.fastq/g' Read_Subsetter.py
sed -i 's/_2.fastq.gz/_R2.fastq/g' Read_Subsetter.py
sed -i 's/fastq.gz/fastq/g' Read_Subsetter.py
```

### Command
```bash
# Organize reads into species-specific folders
mkdir -p Species_A; mv R1_nuclear.fastq R2_nuclear.fastq Species_A/
mkdir -p Species_B; mv R1_nuclear.fastq R2_nuclear.fastq Species_B/

# Subset reads using SISRS read_subsetter.py script
python read_subsetter.py -d Species_A -s 10 -g 500M
python read_subsetter.py -d Species_B -s 10 -g 500M
```

### Output
- **Subsetted FASTQ files** for each species.

---

## **Step 5: Composite Genome Assembly**

### **Step 5.1: Assemble Composite Genome**

### Description
Assemble the composite genome using Ray, ensuring that the output directory starts with "Ray_" as per Ray's requirements.

```bash
# Assemble the composite genome with Ray
mpirun -n 100 Ray -k 31 -detect-sequence-files /projects/SISRS_Project/Composite_Reads -o $output_dir
```

### Output
- **Composite genome assembly** (`composite_genome/Contigs.fasta`)

### **Step 5.2: Prepare Composite Genome**

### Description
Use a custom script to prepare the composite genome, which includes renaming contigs, indexing, and generating site lengths. Ensure directory paths are accurate before running the script.

```bash
./prepare_Composite.sh
```

### Output
- **Prepared composite genome** with contigs renamed and indexed.

---

## **Step 6: Fixed Allele Calling**

### Description
Map each sample to the composite genome to generate sample-specific orthologs. This step only allows for uniquely mapped reads against the composite genome. A strict homozygosity filter is applied based on study requirements. In this study, a minimum of 5 reads for 100% fixed alleles was explicitly specified.

```bash
# Generate SISRS mapping scripts using 64 processors, requiring 5 reads of coverage, and 100% fixed alleles
python Mapping/baseScipts/generate_sisrs_scripts.py 64 5 1

# Navigate to the directory containing the generated scripts
cd /projects/SISRS_Project/test_fastq/Analysis_Dir/Mapping/SISRS_Run/Species_Reads/

# Submit each SLURM script
for script in */*.sh; do
  sbatch $script
done
```

### Output
- **Mapped orthologs** for each sample.

---

## **Step 7: Alignment and Ortholog Generation**

### Description
Use `Output_SISRS.py` to align reads and generate species-specific orthologs. Filter alignments to remove gaps and missing data.

```bash
# Ensure paths are correctly specified in Output_SISRS.py and get_alignment.py
python Output_SISRS.py

# Filter alignments to remove gaps and missing data
python baseScipts/filter_nexus_for_missing.py alignment_singletons.nex

 0
python baseScipts/filter_nexus_for_missing.py alignment_pi.nex 0
```

### Output
- **Alignment files** (`alignment.nex`, `alignment_singletons.nex`, `alignment_pi.nex`, `alignment_pi_singletons.nex`, `alignment_bi.nex`)
- **Filtered alignments** (`alignment_singletons_m0.phylip-relaxed`, `alignment_pi_m0.phylip-relaxed`)

---

## **Step 8: Phylogenetic Analysis for Species Pooling**

### Description
Construct a phylogenetic tree to observe how samples group against the composite genome. Groupings should indicate whether the expected samples are within their clades, ensuring an effective composite genome. Use the parsimony-informative alignment file (`alignment_pi.nex`) for tree generation.

```bash
iqtree -nt AUTO -s alignment_pi_m0.phylip-relaxed -m MFP+MERGE+ASC -bb 1000
```

### Output
- **Phylogenetic tree files** (`.treefile`, `.contree`, `.log`)

---

## **Step 9: Species-Specific SNP Identification**

### Description
Perform all steps for mapping, alignment, and filtering, but this time as pooled samples of a species to generate species-specific fixed alleles. These will be found in the output singletons file.

---

