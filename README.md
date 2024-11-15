# Argo
Argo: species-resolved profiling of **a**ntibiotic **r**esistance **g**enes in complex metagenomes through long-read **o**verlapping

## Introduction
Argo is a long-read-based profiler developed for environmental surveillance of antibiotic resistance genes (ARGs) with species-level resolution. It uses minimap2's base-level alignment with GTDB to obtain raw species assignments and consolidates these assignments on a *read cluster* basis (determined through decomposing the read-overlap graph) by solving a set cover problem. Argo takes *quality-controlled* long reads (either Nanopore or PacBio) as input and returns a table listing predicted ARGs (types and subtypes), their potential hosts, and estimated abundances, expressed as *ARG copies per genome* (cpg), which is equivalent to *ARG copies per cell* (cpc), assuming each cell contains a single genome.


Argo uses [SARG+](https://github.com/xinehc/sarg-curation) as its default database, which augments experimentally validated sequences with RefSeq sequences that share annotation evidence. However, it also accepts customized databases for compatibility with [NDARO](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/) and [CARD](https://card.mcmaster.ca/). See https://github.com/xinehc/argo-supplementary for more details.

## Quick Start
### Installation
Create a new conda environment and install:
```bash
conda create -n argo -c bioconda -c conda-forge argo
conda activate argo
```

### Database setup
Download the database from [Zenodo](https://zenodo.org/records/12576528):
```bash
wget -qN --show-progress https://zenodo.org/records/12576528/files/database.tar.gz
tar -xvf database.tar.gz
```

Index the files:
```bash
## if you encounter memory issue please consider manually lowering cpu_count or simply set cpu_count=1
cpu_count=$(python -c 'import os; print(os.cpu_count())')

diamond makedb --in database/prot.fa --db database/prot --quiet
diamond makedb --in database/sarg.fa --db database/sarg --quiet
ls database/*.*.fa | sort | xargs -P $cpu_count -I {} bash -c '
    filename=${1%.fa*};
    filename=${filename##*/};
    minimap2 -x map-ont -d database/$filename.mmi ${1} 2> /dev/null;
    echo "Indexed <database/$filename.fa>.";' - {}

## remove temporary files to save space
rm -rf database/*.fa
```

### Run Argo
> [!NOTE]
> Argo by default classifies *all reads* that carry at least one ARG into their "most likely" lineages with ties resolved based on the estimated genome copies of species present. Since *plasmid reads* can have multiple hosts in a sample (e.g., [NZ_OW968330.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_OW968330.1)), interpretation requires caution. `--plasmid` forces the splitting of ARGs by their carriers (chromosomes or plasmids), but chimeric reads and uncharacterized plasmids may interfere with the identification.

We provide an example file comprising 10,000 quality-controlled (processed with `Porechop` and `nanoq`) prokaryotic reads (fungal and other reads removed with `minimap2`), randomly selected from the R10.3 mock sample of [Loman Lab Mock Community Experiments](https://lomanlab.github.io/mockcommunity/r10.html).
```bash
wget -qN --show-progress https://zenodo.org/records/12571849/files/example.fa.gz
argo example.fa.gz -d database -o . --plasmid
```

You should see:
```
INFO: Estimating genome copies ...
INFO: ... found 27.375 copies of genomes (bacteria: 27.375; archaea: 0).
INFO: Assigning taxonomy ...
INFO: Reassigning taxonomy ...
INFO: ... found 8 unique species (bacteria: 8; archaea: 0).
INFO: Overlapping ...
INFO: ... median sequence divergence: 0.0519 | initial identity cutoff: 0.9 * 77.03.
INFO: Annotating ARGs ...
INFO: ... candidate HSPs: 11094 | ARG-containing reads: 627.
INFO: Overlapping of ARG-containing reads ...
INFO: ... median sequence divergence of ARG-containing reads: 0.0504 | identity cutoff: 77.40 | low-identity HSPs: 3511.
INFO: Assigning taxonomy ...
INFO: Graph clustering ...
INFO: ... read clusters: 168 | low-subject-cover HSPs: 525 | overlapping HSPs: 6200 | remaining HSPs: 858
INFO: Set covering ...
INFO: ... ARG copies per genome: 28.581 | ARG types: 18 | ARG subtypes: 157 | ARG-carrying species: 8.
INFO: Done.
```

Output file `example.sarg.tsv` lists ARG abundance estimates (cpg) by species:
```text
lineage                     type            subtype    carrier       copy      genome    abundance
...
...Staphylococcus aureus    multidrug       norA       chromosome    4.495     3.875     1.160
...Staphylococcus aureus    multidrug       norB       chromosome    4.380     3.875     1.130
...Staphylococcus aureus    multidrug       norC       chromosome    5.813     3.875     1.500
...Staphylococcus aureus    multidrug       sdrM       chromosome    3.996     3.875     1.031
...Staphylococcus aureus    tetracycline    tet(38)    chromosome    1.749     3.875     0.451
...Staphylococcus aureus    tetracycline    tet(L)     plasmid       14.994    3.875     3.869
...
```

Output file `example.sarg.json` contains detailed annotation information for ARG-containing reads:
```json5
{
    ...
    "dab6192d-f9e3-4eb6-ae8c-3d40e787e072": {
        "remark": "ARG-containing",
        "hit": [
            "multidrug@MFS|norA|WP_001041274.1"
        ],
        "lineage": "Bacteria;Bacillota;Bacilli;Staphylococcales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus",
        "plasmid": false
    },
    ...
    "a9eae701-3767-4a1f-ac00-8945fcde986a": {
        "remark": "ARG-containing",
        "hit": [
            "aminoglycoside|ant(4')-I|WP_001795128.1",
            "tetracycline|tet(L)|WP_012779616.1"
        ],
        "lineage": "Bacteria;Bacillota;Bacilli;Staphylococcales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus",
        "plasmid": true
    },
    ...
}
```

### Run Argo with fixed identity cutoffs
By default, Argo infers the median sequence divergence of a sample using its read overlaps and determines an adaptive identity cutoff for identifying ARGs, calculated as *90 - 2.5 * 100 * median sequence divergence*. This ensures samples generated by different platforms or kits are comparable, despite their differences in read quality. However, if you want to suppress this option and use a fixed cutoff—for instance, an 80% identity cutoff and 80% subject cover—you can run:

```bash
argo *.fa -d database -o . --plasmid -I 80 -S 80
```

A complete list of arguments and their default values is shown below:
```
usage: argo -d DIR -o DIR [-t INT] [-k DIR] [--plasmid] [--skip-clean] [-m INT] [-e FLOAT] [-i FLOAT] [-s FLOAT] [-n INT] [-p FLOAT] [-a INT] [-c FLOAT] [-M INT] [-E FLOAT] [-I FLOAT] [-S FLOAT] [-N INT]
            [-P FLOAT] [-z FLOAT] [-u INT] [-b INT] [-x FLOAT] [-y FLOAT]
            FILE [FILE ...]

Argo: species-resolved profiling of antibiotic resistance genes in complex metagenomes through long-read overlapping

positional arguments:
  FILE                  Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file, gzip optional <*.gz>.

required arguments:
  -d DIR, --db DIR      Unzipped database folder, should contains <prot.fa|sarg.fa>, <nucl.*.fa|sarg.*.fa> and metadata files.
  -o DIR, --output DIR  Output folder.

optional arguments:
  -t INT, --threads INT
                        Number of threads. [10]
  -k DIR, --db-kraken DIR
                        Unzipped kraken2 database for pre-filtering of non-prokaryotic reads. Skip if not given.
  --plasmid             List ARGs carried by plasmids.
  --skip-clean          Skip cleaning, keep all temporary <*.tmp> files.

additional arguments for profiling genomes:
  -m INT                Max. number of target sequences to report (--max-target-seqs/-k in diamond). [25]
  -e FLOAT              Max. expected value to report alignments (--evalue/-e in diamond). [1e-15]
  -i FLOAT              Min. identity in percentage to report alignments (--id in diamond). [0]
  -s FLOAT              Min. subject cover to report alignments (--subject-cover in diamond). [75]
  -n INT                Max. number of secondary alignments to report (-N in minimap2). [2147483647]
  -p FLOAT              Min. secondary-to-primary score ratio to report secondary alignments (-p in minimap2). [0.9]
  -a INT                Terminal condition for EM - max. iterations. [1000]
  -c FLOAT              Terminal condition for EM - epsilon (precision). [1e-10]

additional arguments for profiling antibiotic resistance genes:
  -M INT                Max. number of target sequences to report (--max-target-seqs/-k in diamond). [25]
  -E FLOAT              Max. expected value to report alignments (--evalue/-e in diamond). [1e-5]
  -I FLOAT              Min. identity in percentage to report alignments. If "0" then set 90 - 2.5 * 100 * median sequence divergence. [0]
  -S FLOAT              Min. subject cover of all HPSs within a read cluster to report alignments. [90]
  -N INT                Max. number of secondary alignments to report (-N in minimap2). [2147483647]
  -P FLOAT              Min. secondary-to-primary score ratio to report secondary alignments (-p in minimap2). [0.9]
  -z FLOAT              Min. estimated genome copies of a species to report it ARG copies and abundances. [1]
  -u INT                Max. number of ARG-containing reads per chunk for overlapping. If "0" then use a single chunk. [0]
  -b INT                Graph clustering parameter for MCL - max. iterations. [1000]
  -x FLOAT              Graph clustering parameter for MCL - inflation. [2]
  -y FLOAT              Graph clustering parameter for MCL - expansion. [2]
```


## FAQ
### Does Argo work with long reads of isolates?

Yes, Argo can provide rough estimates of ARG abundances (cpg) for isolates. However, the computational time may be longer for certain pathogenic species (e.g., *Escherichia coli*, *Salmonella enterica*) which typically contain many copies of ARGs on their genomes and are highly redundant in GTDB.

### Why is Argo running slowly for certain samples?

The computational time increases not only with the size of the sample but also with the number of ARG-containing reads and the redundancy of the database. If your sample contains a large proportion of *Escherichia coli* (see above), the computational time is likely to be much longer than usual.

### Does Argo work with assembled contigs?

No, Argo is inherently read-based and does not work with contigs. You may consider using `diamond blastx/p` directly with SARG+/NDARO/CARD for ARG annotation.
