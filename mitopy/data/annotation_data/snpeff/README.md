# Prebuilt database for mitochondrial genome [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/251831106)

## Build process

Download following files from NCBI for [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/251831106) and place to `data/NC_012920` folder:
* complete Genbank record and rename to `genes.gbk`
* complete protein sequence (FASTA) and rename to `protein.fa`

Add record to config file:

```
NC_012920.genome : NC_012920
        NC_012920.chromosomes : chrM
        NC_012920.chrM.codonTable : Vertebrate_Mitochondrial
```

Build database:

```
java -jar snpEff.jar  build -v -genbank NC_012920

```