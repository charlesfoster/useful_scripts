# REROOT MULTIPLE TREES
# Usage: rscript reroot_trees.r outgroup_name
outgroup = commandArgs(trailingOnly=TRUE)

library(ape)

print("Your outgroup is: ")
outgroup

trees <- grep("tre", dir(), value = T)

gene_tree <- vector()

	for(i in 1:length(trees)){
	gene_tree <- read.nexus(trees[i])
	print(paste('I am working on this tree:', trees[i]))
	treeRooted <- root(gene_tree, outgroup, resolve.root=TRUE)
	output_name <- gsub(' ', '_', paste('rooted', trees[i]))
	write.nexus(treeRooted, file = output_name)
}
