# 🧬 Allelic Variant Detection in Parental Lines Using NCBI Sequence Data

This repository provides a complete pipeline for identifying **allelic variants (SNPs/Indels)** present in parental lines using publicly available **NCBI SRA** sequencing data.

---

## 📦 Requirements

Install the following tools before starting:

| Tool         | Use Case                         | Installation |
|--------------|----------------------------------|--------------|
| `sra-tools`  | Download sequencing data         | [Link](https://github.com/ncbi/sra-tools) |
| `bwa`        | Read alignment                   | [Link](http://bio-bwa.sourceforge.net/) |
| `samtools`   | BAM/SAM file handling            | [Link](https://www.htslib.org/) |
| `bcftools`   | Variant calling and comparison   | [Link](https://samtools.github.io/bcftools/) |
| `snpEff`     | (Optional) Variant annotation    | [Link](http://snpeff.sourceforge.net/) |

---

## 🔢 Step 1: Download Sequencing Data from NCBI

1. Go to [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)
2. Search for the parental line and note its **SRA run accession** (e.g., `SRR12345678`)
3. Download the data using `sra-tools`:

```bash
# Download the SRA file
prefetch SRR12345678

# Convert to FASTQ format
fasterq-dump SRR12345678
```

Repeat this for each parental line.

---

## 🌱 Step 2: Download the Reference Genome

1. Go to [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) or [Ensembl Plants](https://plants.ensembl.org/index.html)
2. Search for your species and download the **genomic FASTA (.fna.gz)** file
3. Prepare the reference:

```bash
# Download the genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/XXXXXX/GCA_XXXX_genomic.fna.gz

# Unzip and rename
gunzip GCA_XXXX_genomic.fna.gz
mv GCA_XXXX_genomic.fna reference.fna
```

---

## 🧷 Step 3: Align Reads to the Reference Genome

### 3.1 Index the reference genome

```bash
bwa index reference.fna
```

### 3.2 Align paired-end reads

```bash
bwa mem reference.fna SRR12345678_1.fastq SRR12345678_2.fastq > parent.sam
```

### 3.3 Convert, sort, and index alignment

```bash
samtools view -Sb parent.sam > parent.bam
samtools sort parent.bam -o parent.sorted.bam
samtools index parent.sorted.bam
```

Repeat this for each parental line.

---

## 🧪 Step 4: Call Variants Using `bcftools`

### 4.1 Variant calling

```bash
bcftools mpileup -Ou -f reference.fna parent.sorted.bam | \
bcftools call -mv -Oz -o parent.vcf.gz
```

### 4.2 Index the VCF file

```bash
bcftools index parent.vcf.gz
```

Repeat for each parent.

---

## 🧠 Step 5: (Optional) Annotate Variants with `snpEff`

### 5.1 Run annotation

```bash
snpEff ann reference_db parent.vcf.gz > parent_annotated.vcf
```

### 5.2 To build a custom database (if not available):

```bash
# Place your reference.fna and reference.gff in snpEff/data/reference_db/
# Then build:
snpEff build -gff3 reference_db
```

---

## ⚖️ Step 6: Compare Variants Between Parents

Use `bcftools isec` to identify unique and shared variants:

```bash
bcftools isec parent1.vcf.gz parent2.vcf.gz -p isec_output/
```

This generates:

| File         | Description            |
|--------------|------------------------|
| `0000.vcf`   | Unique to parent1      |
| `0001.vcf`   | Unique to parent2      |
| `0002.vcf`   | Shared between parents |

### (Optional) Variant statistics

```bash
bcftools stats -F reference.fna -s - parent1.vcf.gz parent2.vcf.gz > stats.txt
```

---

## ✅ Output Summary

| File                   | Description                       |
|------------------------|-----------------------------------|
| `parent.sorted.bam`    | Aligned, sorted reads             |
| `parent.vcf.gz`        | Variant calls                     |
| `parent_annotated.vcf` | Annotated variants (optional)     |
| `isec_output/`         | Comparison of parent variants     |

---

## 🧪 Extra Tips

- Always check read quality with `FastQC`
- Automate steps with a Snakemake or bash script
- For structural variation, consider tools like `Delly`, `Manta`, or `LUMPY`

---

## 📘 License

This project is licensed under the [MIT License](LICENSE).

---

## 📬 Contact

Open an issue for questions, contributions, or feedback.
