# Argo
Argo: species-resolved profiling of **a**ntibiotic **r**esistant **g**enes in complex metagenomes through long-read **o**verlapping

## Quick Start
### Installation
```bash
conda crate -n argo -c bioconda -c conda-forge argo
conda activate argo
```

### Database setup
Download:
```bash
wget -qN --show-progress https://zenodo.org/records/12576528/files/database.tar.gz
tar -xvf database.tar.gz
```

Index:
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
We provide an example file comprising 10,000 quality-controlled (processed with `Porechop` and `nanoq`) prokaryotic reads (fungal and other reads removed with `minimap2`), randomly selected from the R10.3 mock sample of [Loman Lab Mock Community Experiments](https://lomanlab.github.io/mockcommunity/r10.html).
```bash
wget -qN --show-progress https://zenodo.org/records/12571849/files/example.fa.gz
argo example.fa.gz -d database -o .
```
