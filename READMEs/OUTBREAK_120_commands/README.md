
##### To turn this data into trees, the gff files generated with prokka are put through roary and iqtree with visualization done with ggtree (in R)
- [roary](https://github.com/sanger-pathogens/Roary)
```
roary -p 48 -f Roary_out -e -n -qc -k kraken_mini_db/minikraken_20141208 Prokka/*/*.gff
```
- [iqtree](http://www.iqtree.org/)
```
iqtree -s Roary_out/core_gene_alignment.aln -t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -pre sample.iqtree -nt 48 -o $control
```
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
