#!/usr/bin/env Rscript

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

# creating the iqtree tree
print("Creating the IQTREE tree")
tree <- read.newick("FULLPATHTOIQTREETREEFILE")
tree_label = get.tree(tree)$tip.label
tree_label_paste <- paste(tree_label,"-UT", sep="")
tree_label_data = data.frame(label=tree_label, label2=substr(tree_label,1,regexpr("-UT",tree_label_paste)-1))
thetree <- ggtree(tree, size=0.25) %<+% tree_label_data + geom_tiplab(color="black", size=1, aes(label=label2), align=TRUE) + geom_treescale(fontsize=1) + labs(title="TITLEOFTREE", caption="DATEOFRUN")
thetree
ggsave("FULLPATHTOTREEIMAGE.pdf", dpi = 900)
ggsave("FULLPATHTOTREEIMAGE_mqc.jpg", dpi = 900)
print("The iqtree tree has been created")

# separating bootstraps for iqtree
bootstrap_values <- thetree$data$label
bootstrap_values[is.na(bootstrap_values)] <- 0
bootstrap_values[bootstrap_values == ""] <- 0
bootstrap_values[which(thetree$data$isTip==TRUE)] <- 0
bootstrap_and_nodes <- data.frame(thetree$data$node, sub("/.*", "", bootstrap_values), sub(".*/", "", bootstrap_values), stringsAsFactors = FALSE)
colnames(bootstrap_and_nodes) <- c("node", "bootstrap","UFboot" )
bootstrap_and_nodes$bootstrap <- as.numeric(bootstrap_and_nodes$bootstrap)
bootstrap_and_nodes$UFboot    <- as.numeric(bootstrap_and_nodes$UFboot)
bootstrapped_tree <- thetree %<+% bootstrap_and_nodes
bootstrap_tree <- bootstrapped_tree + geom_nodepoint(aes(subset = as.numeric (bootstrap) > 70, color = bootstrap), size = 0.5) + scale_color_continuous(low="yellow", high="black") + labs(subtitle="Nodes labeled with bootstrap values") + theme(legend.position="right")
bootstrap_tree
ggsave("FULLPATHTOBOOTSTRAPTREEIMAGE.pdf", dpi = 900)
ggsave("FULLPATHTOBOOTSTRAPTREEIMAGE_mqc.jpg", dpi = 900)
print("Bootstraps have been assigned to iqtree tree")

# creating the gene presence/absence heatmap
print("Reading gene presense/absence table from Roary")
gene_table <- read.table("FULLPATHTOROARYGENETABLE", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
row.names(gene_table) <- gene_table[,1]
gene_table <- gene_table[,-1]
gene_table <- t(gene_table)
print("Gene presense/absence table has been created")

# creating the distance heatmap
print("Calculating pairwise nucleotide distance between tips")
alignment <- read.alignment("FULLPATHTOROARYALIGNMENT", format = "fasta")
alignment_DNAbin <- as.DNAbin(alignment)
alignment_dist <-dist.dna(alignment_DNAbin,variance = TRUE, gamma = TRUE, base.freq = TRUE)
distances <- as.data.frame(as.matrix(alignment_dist))
print("Pairwise nucleotide distances have been calculated")

print("Sorting the nucleotide distances matrix for iqtree tree")
thetree_datasubset <- subset(thetree$data,thetree$data$isTip==TRUE)
thetree_datasubset_sorted <- thetree_datasubset[order(thetree_datasubset$y, decreasing = TRUE),]
sorted_labels <- thetree_datasubset_sorted[,1]
reordered_distances <- distances[sorted_labels]
if (!('control' %in% reordered_distances)) {
  under_max <- max(reordered_distances)
} else {
  under_max <- max(reordered_distances[-grep("control", row.names(reordered_distances)),-grep("control", colnames(reordered_distances))])
}
print("The pairwise nucleotide distances table is complete")

# creating the resistence gene heatmap
print("Creating tables from Abricate resistence results")
resist_table <- read.csv("FULLPATHTOABRICATE_TABLE", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
resist_table_rows <- resist_table[,1]
resist_table_columns <- colnames(resist_table)
resist_table_columns <- resist_table_columns[-1]
resist_table <- as.matrix(resist_table[,-1])
row.names(resist_table) <- resist_table_rows
if ( resist_table_columns[1] == "" ) {
  colnames(resist_table)[1] <- "none"
} else {
  colnames(resist_table) <- resist_table_columns
}
resist_table[is.na(resist_table)] <- 0
resist_table[resist_table == "."] <- 0
if ( ncol(resist_table) == 1 ) {
  resist_table <- resist_table
} else {
  resist_table <- resist_table[,order(colSums(resist_table), decreasing = TRUE)]
}
print("Table of resistence gene presence/absence is complete")

largest_x <- thetree$data$x[order(thetree$data$x, decreasing = TRUE)[1]]
number_of_heatmap <- ncol(distances)
number_of_resist  <- ncol(resist_table)

if ( number_of_heatmap < 25 ) {
distance_size <- 1
} else if ( number_of_heatmap < 50 ) {
distance_size <- 0.75
} else if ( number_of_heatmap < 75 ) {
distance_size <- 0.5
} else {
distance_size <- 0.25
}

if ( number_of_resist < 25 ) {
resist_size <- 1
} else if ( number_of_resist < 50 ) {
resist_size <- 0.75
} else if ( number_of_resist < 75 ) {
resist_size <- 0.5
} else {
resist_size <- 0.25
}

# the tree plus heatmaps
print("Creating heatmap with tree and gene presence/absence table")
gheatmap(thetree + xlim_tree(largest_x * 1.2), gene_table, low="white", high = "black", color = FALSE, colnames=FALSE, offset = largest_x * 0.4) + labs(subtitle="Roary gene presence/absence") + guides(fill=FALSE)
ggsave("FULLPATHTOGENETABLIMAGE.pdf", dpi = 900, width = 4)
ggsave("FULLPATHTOGENETABLIMAGE_mqc.jpg", dpi = 900, width = 4)

print("Creating heatmap with tree and nucleotide pairwise distance matrix")
gheatmap(thetree + xlim_tree(largest_x * 1.2), reordered_distances, color = FALSE, colnames = TRUE, colnames_position = "top", colnames_angle = 90, font.size = distance_size, colnames_offset_y = 0, hjust = 0, offset = largest_x * 0.4) + labs(subtitle="Nucleotide pairwise distance matrix") + scale_fill_gradientn(colors = rev(magma(25)))
ggsave("FULLPATHTODISTANCEIMAGE.pdf", dpi = 900, width = 4)
ggsave("FULLPATHTODISTANCEIMAGE_mqc.jpg", dpi = 900, width = 4)

print("Creating heatmap with tree and nucleotide pairwise distance matrix with limits")
gheatmap(thetree + xlim_tree(largest_x * 1.2), reordered_distances, color = FALSE, colnames = TRUE, colnames_position = "top", colnames_angle = 90, font.size = distance_size, colnames_offset_y = 0, hjust = 0, offset = largest_x * 0.4) + labs(subtitle="Nucleotide pairwise distance matrix") + scale_fill_gradientn(colors = rev(magma(25)), limits = c( 0 , under_max ), na.value="black")
ggsave("FULLPATHTODISTANCEIMAGE_limited.pdf", dpi = 900, width = 4)
ggsave("FULLPATHTODISTANCEIMAGE_limited_mqc.jpg", dpi = 900, width = 4)

print("Creating heatmap with tree and resistence gene presence/absence table")
gheatmap(thetree + xlim_tree(largest_x * 1.2), resist_table, low="white", high = "black", color = FALSE, colnames = TRUE, colnames_position = "top", colnames_angle = 90, font.size = resist_size, colnames_offset_y = 0, hjust = 0, offset = largest_x * 0.4) + labs(subtitle="Resistence gene presence/absence")
ggsave("FULLPATHTORESISTENIMAGE.pdf", dpi = 900, width = 4)
ggsave("FULLPATHTORESISTENIMAGE_mqc.jpg", dpi = 900, width = 4)
print("saved all iqtree images - complete!")
