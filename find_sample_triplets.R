#!/usr/bin/env Rscript
query_triplet <- function(phy, triplet){
  mrca_node <- ape::getMRCA(phy, triplet)
  clade <- treeio::tree_subset(
    tree=phy,
    node=mrca_node,
    levels_back = 0,
    group_node = TRUE,
    group_name = "MATCHED_SUBSET_TREE",
    root_edge = TRUE
  )
  num_tips <- length(clade$tip.label)
  
  if(num_tips == 3){
    trip_clade = TRUE
  } else{
    trip_clade = FALSE
  }
  
  names <- paste(triplet,collapse = "; ")
  result <- data.table::data.table(Sample_Triplet=names,Triplet_Clade=trip_clade,Tips_in_Clade=num_tips, Clade_Support = clade$node.label[1])
  return(result)
}

# read in tree to query
tree <- treeio::read.tree("query.tre") 

# read in samples file, two columns: 
# - first column: header name = tip_name, values must match tip names in tree exactly 
# - second column: header name = sample_identifier, samples in the form of: Sample1_XXXXX, Sample2_XXXX etc., e.g. Sample_1a, Sample_1b, Sample_2a, Sample_2b etc.
# Note: currently restricted to looking for triplets, but the function could be easily modified to look for less or more (e.g., pairs, quads etc.)
sample_table <- data.table::fread("samples.csv",header = TRUE)

samples <- unique(gsub("_.*","_",unlist(sample_table$sample_identifier)))

list_of_sample_triplets <- unique(lapply(1:length(samples), function(x){
  dt <- sample_table %>% 
    filter(grepl(samples[x], sample_identifier,fixed = TRUE))
  triplet <- sort(unique(dt$tip_name))
  return(triplet)
}))

results <- dplyr::bind_rows(lapply(list_of_sample_triplets, query_triplet, phy=tree))
