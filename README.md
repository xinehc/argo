# Argo
Argo: species-resolved profiling of **a**ntibiotic **r**esistant **g**enes in complex metagenomes through long-read **o**verlapping

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
## If you encounter memory issue please consider manually lowering cpu_count or simply set cpu_count=1
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
> Argo by default does not assign taxonomic labels to **plasmid reads** since plasmids can have multiple hosts (for example [NZ_OW968330.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_OW968330.1)). `--plasmid` forces taxonomic classification of these reads by assigning them to their "most likely" lineages (ties resolved based on the estimated genome copies of species within a sample), but interpretation requires caution.

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

Output file `example.sarg.tsv` lists ARG abundance estimates (copy per genome, cpg) by species:
```
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
```
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
