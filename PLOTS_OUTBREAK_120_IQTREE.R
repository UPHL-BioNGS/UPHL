#!/usr/bin/env Rscript

library(ggplot2)
library(ggtree)
library(treeio)
library(seqinr)
library(ape)
library(pegas)
library(ade4)
library(gplots)
library(ggstance)
library(phytools)

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
#bootstrap_and_nodes <- data.frame(thetree$data$node, sub("_.*", "", bootstrap_values), sub(".*_", "", bootstrap_values), stringsAsFactors = FALSE)
colnames(bootstrap_and_nodes) <- c("node", "bootstrap","UFboot" )
bootstrap_and_nodes$bootstrap <- as.numeric(bootstrap_and_nodes$bootstrap)
bootstrap_and_nodes$UFboot    <- as.numeric(bootstrap_and_nodes$UFboot)
bootstrapped_tree <- thetree %<+% bootstrap_and_nodes
bootstrap_tree <- bootstrapped_tree + geom_nodepoint(aes(subset = as.numeric (bootstrap) > 70, color = bootstrap), size = 0.5) + scale_color_continuous(low="yellow", high="black") + labs(subtitle="Nodes labeled with bootstrap values") + theme(legend.position="right")
UFboot_tree <- bootstrapped_tree + geom_nodepoint(aes(subset = as.numeric (UFboot) > 70, color = UFboot), size = 0.5) + scale_color_continuous(low="yellow", high="black") + labs(subtitle="Nodes labeled with UFboot values") + theme(legend.position="right")
bootstrap_tree
ggsave("FULLPATHTOBOOTSTRAPTREEIMAGE.pdf", dpi = 900)
ggsave("FULLPATHTOBOOTSTRAPTREEIMAGE_mqc.jpg", dpi = 900)
UFboot_tree
ggsave("FULLPATHTO_UFBOOT_TREE_IMAGE.pdf", dpi = 900)
ggsave("FULLPATHTO_UFBOOT_TREE_IMAGE_mqc.jpg", dpi = 900)
print("Bootstraps have been assigned to iqtree tree")

# finding distances for iqtree
print("Finding average distance of clade for each node in iqtree tree")
distance_tree <- thetree
df_distance = data.frame()
for(i in 1:tree$Nnode) {
  subtree <- subtrees(tree)[[i]]
 	tree_dist <- cophenetic(subtree)
  mean_clade_distance <- mean(tree_dist)
	tips_for_i <- c(subtrees(tree)[[i]]$tip.label)
	first <- tips_for_i[1]
  last  <- tips_for_i[length(tips_for_i)]
  node_number <- MRCA(thetree, tip = c(first,last))
  df <- data.frame(node_number, mean_clade_distance)
  df_distance <- rbind(df_distance,df)
  }
distance_tree$data$label[which(distance_tree$data$isTip==FALSE)] <- distance_tree$data$node[which(distance_tree$data$isTip==FALSE)]
distance_tree <- distance_tree %<+% df_distance
distance_tree <- distance_tree + geom_nodepoint(aes(color=mean_clade_distance), size=0.5) + scale_color_continuous(low="yellow", high="black") + labs(subtitle="Nodes labeled with average distance in clade") + theme(legend.position="right")
distance_tree
ggsave("FULLPATHTODISTANCETREEIMAGE.pdf", dpi = 900)
ggsave("FULLPATHTODISTANCETREEIMAGE_mqc.jpg", dpi = 900)

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
print("The pairwise nucleotide distances table is complete")

# creating the resistence gene heatmap
print("Creating tables from Abricate resistence results")
resist_table <- read.csv("FULLPATHTOABRICATE_TABLE", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
resist_table_rows <- resist_table[,1]
resist_table <- resist_table[,-1]
resist_table <- resist_table[,-1]
resist_table[is.na(resist_table)] <- 0
resist_table[resist_table == "."] <- 0
resist_table = as.matrix(as.data.frame(lapply(resist_table, as.numeric)))
row.names(resist_table) <- resist_table_rows
resist_table <- resist_table[,order(colSums(resist_table), decreasing = TRUE)]
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
gheatmap(thetree + xlim_tree(largest_x * 1.2), gene_table, low="white", high = "black", color = FALSE, colnames=FALSE, offset = largest_x * 0.4) + labs(subtitle="Roary gene presence/absence")
ggsave("FULLPATHTOGENETABLIMAGE.pdf", dpi = 900, width = 4)
ggsave("FULLPATHTOGENETABLIMAGE_mqc.jpg", dpi = 900, width = 4)

print("Creating heatmap with tree and nucleotide pairwise distance matrix")
gheatmap(thetree + xlim_tree(largest_x * 1.2), reordered_distances, low="red", high = "black", color = FALSE, colnames = TRUE, colnames_position = "top", colnames_angle = 90, font.size = distance_size, colnames_offset_y = 0, hjust = 0, offset = largest_x * 0.4) + labs(subtitle="Nucleotide pairwise distance matrix")
ggsave("FULLPATHTODISTANCEIMAGE.pdf", dpi = 900, width = 4)
ggsave("FULLPATHTODISTANCEIMAGE_mqc.jpg", dpi = 900, width = 4)

print("Creating heatmap with tree and resistence gene presence/absence table")
gheatmap(thetree + xlim_tree(largest_x * 1.2), resist_table, low="white", high = "black", color = FALSE, colnames = TRUE, colnames_position = "top", colnames_angle = 90, font.size = resist_size, colnames_offset_y = 0, hjust = 0, offset = largest_x * 0.4) + labs(subtitle="Resistence gene presence/absence")
ggsave("FULLPATHTORESISTENIMAGE.pdf", dpi = 900, width = 4)
ggsave("FULLPATHTORESISTENIMAGE_mqc.jpg", dpi = 900, width = 4)
print("saved all iqtree images - complete!")
