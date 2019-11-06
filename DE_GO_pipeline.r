#!/usr/bin/env Rscript

################################################################################################
### Purpose: import counts, DE analysis (DESeq2), GO & KEGG enrichment (GOseq) 				 ###
### Most code in here is written by Charles Foster. Some is adapted from the vignettes of    ###
### the various packages used.                                                               ###
################################################################################################

using<-function(...) {
    libs<-unlist(list(...))
    req<-unlist(lapply(libs,require,character.only=TRUE))
    need<-libs[req==FALSE]
    if(length(need)>0){ 
        install.packages(need)
        lapply(need,require,character.only=TRUE)
    }
}

using("argparser")

parser <- arg_parser("DE and GOSeq pipeline. Note that with the exception of '-x', 'optional' arguments aren't actually optional. Please include them.")
parser <- add_argument(parser, c("--samples","--counts_directory","--transcript_to_gene","--gene_lengths","--GO_annotations","--KEGG_annotations","--pathways"), help=c("samples file name", "directory with abundance estimation for each replicate from Salmon","transcript to gene mapping file","gene lengths file","GO annotations for genes","KEGG pathway annotations for genes","KEGG pathway descriptions"), flag=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))
args <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))
if(is.na(args$samples) == T || is.na(args$counts_directory) == T || is.na(args$transcript_to_gene) == T || is.na(args$gene_lengths) == T || is.na(args$GO_annotations) == T || is.na(args$KEGG_annotations) == T || is.na(args$pathways)){
  stop(cat("\nYou have not provided all necessary arguments. \n\nFor help, try: rscript DE_GO_pipeline.r -h\n"), call.=FALSE)
} 

using("tximport")
using("DESeq2")
using("goseq")
using("GO.db")
using("edgeR")
using("readr")
using("qdapTools")
using("dplyr")
using("RColorBrewer")
using("pheatmap")
using("ggplot2")

###### Step 1: Import the count data into a tximport object ######
cat("\n*** Importing the data ***\n")
dir <- getwd()
cat("\nWorking in:", dir,"\n")
samples <- read.table(args$samples,header=T,stringsAsFactors=F)
rownames(samples) <- samples$Replicate

num_g <- length(unique(samples[,1]))

expand.grid.unique <- function(x, y, include.equals=FALSE)
{
    x <- unique(x)
    y <- unique(y)
    g <- function(i)
    {
        z <- setdiff(y, x[seq_len(i-include.equals)])
        if(length(z)) cbind(x[i], z, deparse.level=0)
    }
    do.call(rbind, lapply(seq_along(x), g))
}

cont.mat <- expand.grid.unique(samples[,1],samples[,1])
contrasts <- paste(cont.mat[,1],cont.mat[,2],sep="_vs_")

cat("\nYou have ",num_g," groups.\n")
cat("\nYour groups are ",unique(samples[,1]),"\n")
cat("\nYour contrasts are ",contrasts,"\n")

reslist = list.files(pattern = "*UP.subset")

if(length(reslist) == 0){
	files <- file.path(args$counts_directory, samples$Replicate,"quant.sf")
	names(files) <- basename(file.path(args$counts_directory, samples$Replicate))
	if(all(file.exists(files)) == TRUE){
	cat("All files exist, and will be read in properly.\n")
	}
	tx2gene <- read.delim(file.path(args$transcript_to_gene))
	cat("\nYou have chosen to use the standard DESeq2 pipeline\n\n")
	txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
	
	group <- as.factor(sort(samples$Group))
	
	###### Step 2: Run DESeq2 to get DE genes ######
	cat("\n*** Running DESeq2: Finding DE genes ***\n\n")
	
	dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~Group)
	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep,]
	
	dds <- DESeq(dds)
	writeLines(rownames(counts(dds)),con="filtered_genes_list.txt")
	
	vsd <- vst(dds, blind=FALSE)
	vsd_counts <- assay(vsd)
	write.table(vsd_counts, file=paste0("vsd_normalised_counts.txt"), sep='\t', quote=FALSE)
	sampleDists <- dist(t(assay(vsd)))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(vsd$Group)
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	pdf("sample_correlation.pdf")
	pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
	dev.off()
	plotPCA(vsd, intgroup=c("Group"))
	ggsave("PCA_plot.pdf")
	
	for(i in 1:length(cont.mat[,1])){
		contrast <- paste(cont.mat[i,1],cont.mat[i,2],sep="_vs_")
		cat("\nWorking with: ",contrast,"\n")
		result <- results(dds, altHypothesis = "greaterAbs", contrast=c("Group",cont.mat[i,1],cont.mat[i,2]), lfcThreshold=1,alpha=0.001)
		baseMeanA_result <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$Group == cont.mat[i,1]])
		baseMeanB_result <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$Group == cont.mat[i,2]])
		result = cbind(baseMeanA_result,baseMeanB_result,as.data.frame(result))
		result = cbind(sampleA=cont.mat[i,1], sampleB=cont.mat[i,2], as.data.frame(result))
		result  <- result[order(result$pvalue),]
		result_sig <- subset(result, padj < 0.001)
		result_UP1 <- result_sig[c(which(result_sig$log2FoldChange >=1)),]
		result_UP2 <- result_sig[c(which(result_sig$log2FoldChange <=1)),]
	
		### Print summary
		cat("\nYour results are... *drumroll*\n")
		cat(contrast,nrow(result_sig),"DE; ",nrow(result_UP1),"upregulated in ",cont.mat[i,1],";",nrow(result_UP2),"upregulated in ",cont.mat[i,2],".")
	
		### Get count data
		ind <- c(grep(cont.mat[i,1],colnames(counts(dds))),grep(cont.mat[i,2],colnames(counts(dds))))
		cts_result = counts(dds)[,ind,drop=F]
		
		### Write DE results
		write.table(result,  file=paste0(cont.mat[i,1],"_vs_",cont.mat[i,2],".DESeq2",".all_genes"), sep='\t', quote=FALSE)
		write.table(result_sig,  file=paste0(cont.mat[i,1],"_vs_",cont.mat[i,2],".DESeq2",".DE_results"), sep='\t', quote=FALSE) 
		write.table(result_UP1, file=paste0(cont.mat[i,1],"_vs_",cont.mat[i,2],".DESeq2",".DE_results.P0.001_C0.",cont.mat[i,1],"-UP.subset"), sep='\t', quote=FALSE)
		write.table(result_UP2, file=paste0(cont.mat[i,1],"_vs_",cont.mat[i,2],".DESeq2",".DE_results.P0.001_C0.",cont.mat[i,2],"-UP.subset"), sep='\t', quote=FALSE)
		
		write.table(cts_result, file=paste0(cont.mat[i,1],"_vs_",cont.mat[i,2],".DESeq2",".count_matrix"), sep='\t', quote=FALSE)
		}
	
	cat("\n\n*** Finished with DESeq2 ***\n\n")
		### Prepare GOSeq files
	### Note: a lot of code in the next part is taken/adapted from the Trinity GOseq script ###
	cat("\n*** Running GOseq: Finding enriched/depleted GO categories/KEGG pathways ***\n\n")
	gene_lengths = read.table(args$gene_lengths, header=T, row.names=1, com='')
	gene_lengths = as.matrix(gene_lengths[,1,drop=F])
	sample_set_gene_ids = readLines("filtered_genes_list.txt")
	sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]
	
	cat("\n*** Preparing the GO annotations. Be patient. ***\n\n")
	GO_info = read.table(args$GO_annotations, header=F, row.names=1,stringsAsFactors=F)
	GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
	names(GO_info_listed) = rownames(GO_info)
	
	get_GO_term_descr =  function(x) {
	    d = 'none';
	    go_info = GOTERM[[x]];
	    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
	    return(d);
	}
	
	GO_to_gene_list = list()
	for (gene_id in intersect(names(GO_info_listed), sample_set_gene_ids)) {
	    go_list = GO_info_listed[[gene_id]]
	    for (go_id in go_list) {
	        GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
	    }
	}
	
	all_gene_ids = rep(0,length(names(sample_set_gene_lengths)))
	names(all_gene_ids) <- names(sample_set_gene_lengths)
	
	cat("\n*** Preparing the KEGG pathways. Be patient. ***\n\n")
	#### parse KEGG pathways
	##### Note: I am now using KEGG pathways sourced from KO numbers, which were originally 
	##### identified by Trinotate. To source these KEGG pathways I used a separate script with KEGGlink
	##### THIS IS NOT STRAIGHT FORWARD! If too hard, comment out everything to do with KEGG.
	
	tmp0 <- read.delim(args$KEGG_annotations,header=F,sep="\t",stringsAsFactors=F)
	tmp1 <- scan(args$KEGG_annotations, what="", sep="\n")
	tmp2 <- strsplit(tmp1, "[[:space:]]+")
	names(tmp2) <- sapply(tmp2, `[[`, 1)
	kegg_anno <- lapply(tmp2, `[`, -1)
	kegg_pathway_names <- read.delim(args$pathways,colClasses=c("numeric",rep("character",3)))
	kegg_split<-strsplit(tmp0$V2," ")
	names(kegg_split)<-tmp0$V1
	kegg_annotation <- lapply(tmp2, `[`, -1)
	kegg_anno_gene <- list2df(kegg_annotation,col1 = "KEGG_ID",col2="GENE")
	
	
	filelist = list.files(pattern = "*UP.subset")
	
	for(file in filelist){
		cat("\nWorking with",basename(file),"\n")
		old_name <- file
		go_enrich_filename <- gsub("subset","enriched.GO.txt",old_name)
		go_depleted_filename <- gsub("subset","depleted.GO.txt",old_name)
		KEGG_enrich_filename <- gsub("subset","enriched.KEGG.txt",old_name)
		KEGG_depleted_filename <- gsub("subset","depleted.KEGG.txt",old_name)
		contrast <- read.table(file,header=T,row.names=1)
		names_contrast <- rownames(contrast)
		ind <- which(names(all_gene_ids) %in% names_contrast)
		contrast.goseq <- all_gene_ids
		contrast.goseq[ind] <- 1
		gene_ids_in_feature_cat = names(contrast.goseq[ind])
		if(all(names(contrast.goseq) != names(sample_set_gene_lengths))){
	    stop("Error, your lengths are jumbled up. Please get this bit right, moron.")}
		pwf = nullp(contrast.goseq, bias.data=sample_set_gene_lengths)
		rownames(pwf) = sample_set_gene_ids
		cat("\nFirstly, I'll search for enriched/depleted GO terms\n")
		contrast.result.GO = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=F)
		contrast.enriched.GO <- contrast.result.GO
		contrast.depleted.GO <- contrast.result.GO
		contrast.enriched.GO$over_represented_FDR = p.adjust(contrast.enriched.GO$over_represented_pvalue, method="BH")
		result_table = contrast.enriched.GO[contrast.enriched.GO$over_represented_FDR<=0.05,]
		descr = unlist(lapply(result_table$category, get_GO_term_descr))
		result_table$go_term = descr;
		result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
		            gene_list = GO_to_gene_list[[x]]
		            gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
		            paste(gene_list, collapse=', ');
		     }) )
		enriched_count <- length(result_table$category)
		if((length(enriched_count) <1) == T){
			cat("\nThere are no enriched GO terms. I'm really sorry.\n")} else {
			cat("\n~~~~~ There are",enriched_count,"enriched GO terms.")
			write.table(result_table[order(result_table$over_represented_FDR),],file=go_enrich_filename, sep='\t', quote=F, row.names=F)
		}
		contrast.depleted.GO$under_represented_FDR = p.adjust(contrast.depleted.GO$under_represented_pvalue, method="BH")
		result_table = contrast.depleted.GO[contrast.depleted.GO$under_represented_FDR<=0.05,]
		descr = unlist(lapply(result_table$category, get_GO_term_descr))
		result_table$go_term = descr;
		result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
		            gene_list = GO_to_gene_list[[x]]
		            gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
		            paste(gene_list, collapse=', ');
		     }) )
		depleted_count <- length(result_table$category)
		if((length(depleted_count) <1) == T){
		cat("\nThere are no depleted GO terms. I'm really sorry.\n")} else {
		cat(" There are",depleted_count,"depleted GO terms. ~~~~~\n")
		write.table(result_table[order(result_table$under_represented_FDR),],file=go_depleted_filename, sep='\t', quote=F, row.names=F)}
	
		cat("\nSecondly, I'll search for enriched/depleted KEGG pathways\n")
		contrast.result.KEGG = goseq(pwf,gene2cat=kegg_anno, use_genes_without_cat=F)
		contrast.enriched.KEGG <- contrast.result.KEGG
		contrast.depleted.KEGG <- contrast.result.KEGG
	
		contrast.enriched.KEGG$over_represented_FDR = p.adjust(contrast.enriched.KEGG$over_represented_pvalue, method="BH")
		result_table = contrast.enriched.KEGG[contrast.enriched.KEGG$over_represented_FDR<=0.05,]
		result_table$description <- rep("",length(result_table$category))
		for(i in 1:length(result_table$description)){
			result_table$description[i] <- kegg_pathway_names$description[kegg_pathway_names$reduced_path_id %in% result_table$category[i]]}
			result_table <- result_table[,c(1,7,2,3,4,5,6)]
			res_vec <- result_table$category
			if((length(res_vec) <1) == T){
			cat("\nThere are no enriched KEGG pathways. I'm really sorry.\n")} else {
			result_table$gene_ids <- rep("",length(res_vec))
			for(i in 1:length(res_vec)){
				result_table$gene_ids[i] <- paste(kegg_anno_gene[kegg_anno_gene$KEGG_ID %in% res_vec[i],][,2],collapse=", ")
			}
			enriched_count <- length(result_table$category)
			cat("\n~~~~~ There are",enriched_count,"enriched KEGG pathways.")
			write.table(result_table[order(result_table$over_represented_FDR),],file=KEGG_enrich_filename, sep='\t', quote=F, row.names=F)
		}
		
		contrast.depleted.KEGG$under_represented_FDR = p.adjust(contrast.depleted.KEGG$under_represented_pvalue, method="BH")
		result_table = contrast.depleted.KEGG[contrast.depleted.KEGG$under_represented_FDR<=0.05,]
		result_table$description <- rep("",length(result_table$category))
		for(i in 1:length(result_table$description)){
			result_table$description[i] <- kegg_pathway_names$description[kegg_pathway_names$reduced_path_id %in% result_table$category[i]]}
			result_table <- result_table[,c(1,7,2,3,4,5,6)]
			res_vec <- result_table$category
			if((length(res_vec) <1) == T){
			cat("\nThere are no depleted KEGG pathways. I'm really sorry.\n")} else {
			result_table$gene_ids <- rep("",length(res_vec))
			for(i in 1:length(res_vec)){
				result_table$gene_ids[i] <- paste(kegg_anno_gene[kegg_anno_gene$KEGG_ID %in% res_vec[i],][,2],collapse=", ")
			}
			enriched_count <- length(result_table$category)
			cat(" There are",enriched_count,"depleted KEGG pathways. ~~~~~\n")
			write.table(result_table[order(result_table$under_represented_FDR),],file=KEGG_depleted_filename, sep='\t', quote=F, row.names=F)

		}
	}
				cat("\n\n*** Your analysis is finished ***\n\n")
} else{
	cat("\n\n*** You have already run DESeq2. Skipping straight to GOSeq ***\n\n")
	### Prepare GOSeq files
	### Note: a lot of code in the next part is taken/adapted from the Trinity GOseq script ###
	cat("\n*** Running GOseq: Finding enriched/depleted GO categories/KEGG pathways ***\n\n")
	gene_lengths = read.table(args$gene_lengths, header=T, row.names=1, com='')
	gene_lengths = as.matrix(gene_lengths[,1,drop=F])
	sample_set_gene_ids = readLines("filtered_genes_list.txt")
	sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]
	
	cat("\n*** Preparing the GO annotations. Be patient. ***\n\n")
	GO_info = read.table(args$GO_annotations, header=F, row.names=1,stringsAsFactors=F)
	GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
	names(GO_info_listed) = rownames(GO_info)
	
	get_GO_term_descr =  function(x) {
	    d = 'none';
	    go_info = GOTERM[[x]];
	    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
	    return(d);
	}
	
	GO_to_gene_list = list()
	for (gene_id in intersect(names(GO_info_listed), sample_set_gene_ids)) {
	    go_list = GO_info_listed[[gene_id]]
	    for (go_id in go_list) {
	        GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
	    }
	}
	
	all_gene_ids = rep(0,length(names(sample_set_gene_lengths)))
	names(all_gene_ids) <- names(sample_set_gene_lengths)
	
	cat("\n*** Preparing the KEGG pathways. Be patient. ***\n\n")
	#### parse KEGG pathways
	##### Note: I am now using KEGG pathways sourced from KO numbers, which were originally 
	##### identified by Trinotate. To source these KEGG pathways I used a separate script with KEGGlink
	##### THIS IS NOT STRAIGHT FORWARD! If too hard, comment out everything to do with KEGG.
	
	tmp0 <- read.delim(args$KEGG_annotations,header=F,sep="\t",stringsAsFactors=F)
	tmp1 <- scan(args$KEGG_annotations, what="", sep="\n")
	tmp2 <- strsplit(tmp1, "[[:space:]]+")
	names(tmp2) <- sapply(tmp2, `[[`, 1)
	kegg_anno <- lapply(tmp2, `[`, -1)
	kegg_pathway_names <- read.delim(args$pathways,colClasses=c("numeric",rep("character",3)))
	kegg_split<-strsplit(tmp0$V2," ")
	names(kegg_split)<-tmp0$V1
	kegg_annotation <- lapply(tmp2, `[`, -1)
	kegg_anno_gene <- list2df(kegg_annotation,col1 = "KEGG_ID",col2="GENE")
	
	
	filelist = list.files(pattern = "*UP.subset")
	
	for(file in filelist){
		cat("\nWorking with",basename(file),"\n")
		old_name <- file
		go_enrich_filename <- gsub("subset","enriched.GO.txt",old_name)
		go_depleted_filename <- gsub("subset","depleted.GO.txt",old_name)
		KEGG_enrich_filename <- gsub("subset","enriched.KEGG.txt",old_name)
		KEGG_depleted_filename <- gsub("subset","depleted.KEGG.txt",old_name)
		contrast <- read.table(file,header=T,row.names=1)
		names_contrast <- rownames(contrast)
		ind <- which(names(all_gene_ids) %in% names_contrast)
		contrast.goseq <- all_gene_ids
		contrast.goseq[ind] <- 1
		gene_ids_in_feature_cat = names(contrast.goseq[ind])
		if(all(names(contrast.goseq) != names(sample_set_gene_lengths))){
	    stop("Error, your lengths are jumbled up. Please get this bit right, moron.")}
		pwf = nullp(contrast.goseq, bias.data=sample_set_gene_lengths)
		rownames(pwf) = sample_set_gene_ids
		cat("\nFirstly, I'll search for enriched/depleted GO terms\n")
		contrast.result.GO = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=F)
		contrast.enriched.GO <- contrast.result.GO
		contrast.depleted.GO <- contrast.result.GO
		contrast.enriched.GO$over_represented_FDR = p.adjust(contrast.enriched.GO$over_represented_pvalue, method="BH")
		result_table = contrast.enriched.GO[contrast.enriched.GO$over_represented_FDR<=0.05,]
		descr = unlist(lapply(result_table$category, get_GO_term_descr))
		result_table$go_term = descr;
		result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
		            gene_list = GO_to_gene_list[[x]]
		            gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
		            paste(gene_list, collapse=', ');
		     }) )
		enriched_count <- length(result_table$category)
		if((length(enriched_count) <1) == T){
			cat("\nThere are no enriched GO terms. I'm really sorry.\n")} else {
			cat("\n~~~~~ There are",enriched_count,"enriched GO terms.")
			write.table(result_table[order(result_table$over_represented_FDR),],file=go_enrich_filename, sep='\t', quote=F, row.names=F)
		}
		contrast.depleted.GO$under_represented_FDR = p.adjust(contrast.depleted.GO$under_represented_pvalue, method="BH")
		result_table = contrast.depleted.GO[contrast.depleted.GO$under_represented_FDR<=0.05,]
		descr = unlist(lapply(result_table$category, get_GO_term_descr))
		result_table$go_term = descr;
		result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
		            gene_list = GO_to_gene_list[[x]]
		            gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
		            paste(gene_list, collapse=', ');
		     }) )
		depleted_count <- length(result_table$category)
		if((length(depleted_count) <1) == T){
		cat("\nThere are no depleted GO terms. I'm really sorry.\n")} else {
		cat(" There are",depleted_count,"depleted GO terms. ~~~~~\n")
		write.table(result_table[order(result_table$under_represented_FDR),],file=go_depleted_filename, sep='\t', quote=F, row.names=F)}
	
		cat("\nSecondly, I'll search for enriched/depleted KEGG pathways\n")
		contrast.result.KEGG = goseq(pwf,gene2cat=kegg_anno, use_genes_without_cat=F)
		contrast.enriched.KEGG <- contrast.result.KEGG
		contrast.depleted.KEGG <- contrast.result.KEGG
	
		contrast.enriched.KEGG$over_represented_FDR = p.adjust(contrast.enriched.KEGG$over_represented_pvalue, method="BH")
		result_table = contrast.enriched.KEGG[contrast.enriched.KEGG$over_represented_FDR<=0.05,]
		result_table$description <- rep("",length(result_table$category))
		for(i in 1:length(result_table$description)){
			result_table$description[i] <- kegg_pathway_names$description[kegg_pathway_names$reduced_path_id %in% result_table$category[i]]}
			result_table <- result_table[,c(1,7,2,3,4,5,6)]
			res_vec <- result_table$category
			if((length(res_vec) <1) == T){
			cat("There are no enriched KEGG pathways. I'm really sorry.\n")} else {
			result_table$gene_ids <- rep("",length(res_vec))
			for(i in 1:length(res_vec)){
				result_table$gene_ids[i] <- paste(kegg_anno_gene[kegg_anno_gene$KEGG_ID %in% res_vec[i],][,2],collapse=", ")
			}
			enriched_count <- length(result_table$category)
			cat("\n~~~~~ There are",enriched_count,"enriched KEGG pathways.")
			write.table(result_table[order(result_table$over_represented_FDR),],file=KEGG_enrich_filename, sep='\t', quote=F, row.names=F)
		}
		
		contrast.depleted.KEGG$under_represented_FDR = p.adjust(contrast.depleted.KEGG$under_represented_pvalue, method="BH")
		result_table = contrast.depleted.KEGG[contrast.depleted.KEGG$under_represented_FDR<=0.05,]
		result_table$description <- rep("",length(result_table$category))
		for(i in 1:length(result_table$description)){
			result_table$description[i] <- kegg_pathway_names$description[kegg_pathway_names$reduced_path_id %in% result_table$category[i]]}
			result_table <- result_table[,c(1,7,2,3,4,5,6)]
			res_vec <- result_table$category
			if((length(res_vec) <1) == T){
			cat("There are no depleted KEGG pathways. I'm really sorry.\n")} else {
			result_table$gene_ids <- rep("",length(res_vec))
			for(i in 1:length(res_vec)){
				result_table$gene_ids[i] <- paste(kegg_anno_gene[kegg_anno_gene$KEGG_ID %in% res_vec[i],][,2],collapse=", ")
			}
			enriched_count <- length(result_table$category)
			cat(" There are",enriched_count,"depleted KEGG pathways. ~~~~~\n")
			write.table(result_table[order(result_table$under_represented_FDR),],file=KEGG_depleted_filename, sep='\t', quote=F, row.names=F)
		}
	}
			cat("\n\n*** Your analysis is finished ***\n\n")
}