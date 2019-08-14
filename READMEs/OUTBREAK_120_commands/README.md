# To detect outbreaks, gff files are organized by species into folders

Mash is used to determine the species of each sample in the UPHL Reference-Free Workflow. This result is used to organize gff into file structures like the following:
`ecoli/O157/H7/*gff`

Then, the following roary and iqtree commands are performed on those directories:

- [roary](https://github.com/sanger-pathogens/Roary)
```
roary -p 48 -f Roary_out -e -n -qc -k kraken_mini_db/minikraken_20141208 ecoli/O157/H7/*gff
```
- [iqtree](http://www.iqtree.org/)
```
iqtree -s Roary_out/core_gene_alignment.aln -t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -pre sample.iqtree -nt 48 -o $control
```

A tree is then generated with ggtree with additional information organized by ape with the following R libraries
- [ape](https://cran.r-project.org/web/packages/ape/index.html) & [ggtree](http://bioconductor.org/packages/release/bioc/html/ggtree.html) (For more information, see [PLOTS_IQTREE.R](outbreak_120_scripts/PLOTS_IQTREE.R))
```
library(ggplot2)
library(ggtree)
library(treeio)
library(seqinr)
library(ape)
library(ade4)
library(gplots)
library(ggstance)
library(phytools)
library(viridis)
```
