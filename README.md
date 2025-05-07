# Allelic-Variant-Detection
Allelic Variant Detection in Parental Lines Using NCBI Sequence Data
Allelic Variant Detection in Parental Lines Using NCBI Sequence Data
This pipeline guides you through identifying allelic variants (SNPs/Indels) present in parental lines using publicly available sequence data from the NCBI Sequence Read Archive (SRA).

📦 Requirements
Make sure the following tools are installed and accessible in your environment:

sra-tools (for downloading FASTQ files)

bwa (for alignment)

samtools (for processing BAM files)

bcftools (for variant calling and comparison)

snpEff (optional, for variant annotation)

🔍 Step 1: Download Sequence Data from NCBI SRA
First, find the parental lines on NCBI SRA and copy the SRA accession number (e.g., SRR12345678).

Download FASTQ files:

bash
Copy
Edit
# Download SRA file
prefetch SRR12345678

# Convert SRA to FASTQ
fasterq-dump SRR12345678
Repeat for all parental samples.

📂 Step 2: Download the Reference Genome
Get a high-quality reference genome (FASTA format) from:

NCBI Assembly

Ensembl Plants

Example:

bash
Copy
Edit
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/XXXXXX/.../GCA_XXXX_genomic.fna.gz
gunzip GCA_XXXX_genomic.fna.gz
🧬 Step 3: Align Reads to the Reference Genome
Index the reference:

bash
Copy
Edit
bwa index reference.fna
Align reads for each parent:

bash
Copy
Edit
# Example for paired-end reads
bwa mem reference.fna parent_R1.fastq parent_R2.fastq > parent.sam
Convert SAM to sorted BAM:

bash
Copy
Edit
samtools view -Sb parent.sam | samtools sort -o parent.bam
samtools index parent.bam
🧪 Step 4: Call Variants with bcftools
bash
Copy
Edit
# Create BCF and call variants
bcftools mpileup -Ou -f reference.fna parent.bam | \
bcftools call -mv -Oz -o parent.vcf.gz

# Index the VCF
bcftools index parent.vcf.gz
Repeat this for each parent.

🧠 Step 5 (Optional): Annotate Variants with snpEff
If you want to understand the functional impact of the variants:

bash
Copy
Edit
# Run snpEff annotation
snpEff ann your_genome_db parent.vcf.gz > annotated_parent.vcf
Replace your_genome_db with the appropriate SnpEff genome name. You may need to build your own database if it's not available.

⚖️ Step 6: Compare Variants Between Parents
Use bcftools isec or bcftools compare to find shared and unique alleles:

bash
Copy
Edit
# Compare VCFs from Parent1 and Parent2
bcftools isec parent1.vcf.gz parent2.vcf.gz -p isec_output/
This will generate:

Shared SNPs/Indels

Unique variants in each parent

Visual Venn-style comparison files

📁 Output Files
File/Folder	Description
parent.bam	Sorted alignment file
parent.vcf.gz	Variant calls (compressed)
annotated.vcf	Annotated variants (optional)
isec_output/	Directory with comparison results

📝 Notes
Ensure read quality by optionally running FastQC before alignment.

Use multiqc for aggregated quality reports.

Consider using GATK or DeepVariant for more advanced variant calling.

🧬 Example
If you have the following data:

Parents: SRR12345678 (Parent1), SRR87654321 (Parent2)

Reference: brassica_ref.fna

The steps would be:

bash
Copy
Edit
# Download and convert
fasterq-dump SRR12345678
fasterq-dump SRR87654321

# Align
bwa mem brassica_ref.fna SRR12345678_1.fastq SRR12345678_2.fastq | samtools sort -o parent1.bam
bwa mem brassica_ref.fna SRR87654321_1.fastq SRR87654321_2.fastq | samtools sort -o parent2.bam

# Index BAM
samtools index parent1.bam
samtools index parent2.bam

# Call variants
bcftools mpileup -Ou -f brassica_ref.fna parent1.bam | bcftools call -mv -Oz -o parent1.vcf.gz
bcftools mpileup -Ou -f brassica_ref.fna parent2.bam | bcftools call -mv -Oz -o parent2.vcf.gz

# Compare
bcftools isec parent1.vcf.gz parent2.vcf.gz -p isec_output/
📬 Contact
If you find this helpful or need support, feel free to open an issue or fork the repository.

📘 License
This pipeline is released under the MIT License.
