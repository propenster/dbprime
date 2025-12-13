# dbprime - a toolkit for finding InDels and developing Molecular markers
dbprime is a Python/C++ toolkit for primer design, and selection of molecular markers (SNPs / InDels) with a focus on maize breeding and comparative genomics.


* Fast (built on minimap2, samtools, Primer3)
* Flexible (pairwise or multi-sample analysis)
* Breeder-friendly (outputs primer-ready marker candidates)
* Reproducible (Conda/Bioconda environment provided)


---

What dbprime Does

Given one reference genome/assembly and one or more query assemblies or FASTA/FASTQ files, dbprime will:

. Align query sequences to the reference (via minimap)
. Detect SNPs and InDels (with emphasis on InDels)
. Extract flanking sequences around each variant
. Design PCR primers flanking each candidate marker (via Primer)
. Output structured results for downstream genotyping or breeding pipelines

Typical applications:

* Maize cultivar comparison
* Marker-assisted selection (MAS)
* Population genotyping
* Comparative genomics

---

Repository Structure

```
dbprime/
├── dbprime.py               Main CLI entry point
├── primer_design.py         Core logic (alignment, indel calling, primer design)
├── examples/
│   └── dbprime/             Example input/output data
├── env.yaml                 Bioconda environment definition
├── README.md                This file
```

---

Installation

You have two options: recommended Bioconda/Conda installation, or with pip.

---

Option : Recommended (Bioconda / Conda)

This ensures all compiled bioinformatics tools work correctly.

. Create environment

```bash
conda env create -f env.yaml
conda activate mindel
```

. Environment definition (`env.yaml`)

```yaml
name: mindel
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=.

   #alignment tools
  - minimap
  - samtools
  - seqtk

   #python libraries
  - biopython
  - pysam
  - pandas
  - primer-py

  - pip:
      - tqdm
      - rich
```

---

Option : pip (Not recommended, but possible)

You must install system binaries manually first:

* minimap2
* samtools
* primer3

Then:

```bash
pip install biopython pysam pandas primer-py tqdm rich
```

This may fail on some systems due to compiled dependencies.

---

Usage

The main entry point is `dbprime.py`.

Pairwise comparison (two samples)

```bash
python dbprime.py pair \
  -q query.fasta \
  -r reference.fasta \
  -o results_pair/ \
  -l  \
  -f  \
  -s   \
  -p  \
  -t 
```

Multi-sample comparison (population mode)

```bash
python dbprime.py multi \
  -q queries_dir/ \
  -r reference.fasta \
  -o results_population/ \
  -l  \
  -f  \
  -s   \
  -p  \
  -t 
```

---

Command-line Arguments

Common arguments

| Flag                       | Description                         | Default  |
| -------------------------- | ----------------------------------- | -------- |
| `-q, --query / --queries`  | Query FASTA/FASTQ file or directory | required |
| `-r, --reference`          | Reference FASTA file                | required |
| `-o, --output`             | Output directory                    | required |
| `-l, --min_indel_length`   | Minimum indel length                |          |
| `-f, --flank`              | Flanking sequence length            |          |
| `-s, --product_size_range` | PCR product size range              | –        |
| `-p, --num_primers`        | Primer pairs per marker             |          |
| `-t, --threads`            | Threads for alignment               |          |

---

Outputs

Each run produces:

* `indel_candidates.tsv` – tabular list of markers
* `indel_candidates.json` – full structured metadata
* Primer sequences (LEFT / RIGHT)
* Coordinates relative to reference genome

These outputs are directly usable for:

* PCR validation
* Marker-assisted breeding
* Database ingestion

---

Maize Breeding Context

dbprime was designed with plant breeding workflows in mind:

* Works with draft or polished assemblies
* Handles cultivar-to-cultivar comparisons
* Supports population-level marker discovery
* Produces primer-ready outputs for wet-lab validation

---
