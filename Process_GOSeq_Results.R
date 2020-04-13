#!/usr/bin/env Rscript

################################################################################################
### Timing: use this script after DE_GO_pipeline.r                          				 ###
###         																				 ###
### Usage: rscript Process_GOSeq_Results.R <gene_mapping_file> <go_annotations_file 		 ###
###         																				 ###
### Purpose: processes GOSeq results to add gene symbols to the enriched/depleted GO terms,  ###
###          writes the resulting tables to file, then plots the results for                 ###
###          enriched/depleted biological process terms as an 'emapplot' and 'dotplot'.      ###
###         																				 ###
### Authorship: Essentially all code in here is written by Charles Foster.                   ###
################################################################################################

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(qdapTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

dir <- getwd()
subDir_enriched_plots <- paste0(dir,"/Enriched/Plots/")
subDir_enriched_tables <- paste0(dir,"/Enriched/Tables/")
subDir_depleted_plots <- paste0(dir,"/Depleted/Plots/")
subDir_depleted_tables <- paste0(dir,"/Depleted/Tables/")
subfolder_names <- c("Enriched/Tables", "Enriched/Plots", "Depleted/Tables", "Depleted/Plots")
cat("\nFirst I'll create results directories...\n")
sapply(subfolder_names, function(x){dir.create(file.path(dir, x),recursive=T)})

mapping <- read.delim(args[1],stringsAsFactors=F,header=F)
mapping <- mapping[which(!is.na(mapping$V2)),]
mapping <- unique(mapping)
mapping <- mapping[!duplicated(mapping$V1),]
mapping <- mapping[order(mapping$V2),]
mapping$V2 <- toupper(mapping$V2)
mapping$V2 <- make.unique(mapping$V2,sep=".")

### Prepare the background: GO
goterms_term <- as.data.frame(Term(GOTERM))
goterms_term$GO <- rownames(goterms_term)
goterms_ont <- as.data.frame(Ontology(GOTERM))
goterms_ont$GO <- rownames(goterms_ont)
goterms<-left_join(goterms_term,goterms_ont,by = "GO")
colnames(goterms) <- c("NAMES","GO","ONTOLOGY")
tmp1 <- scan(args[2], what="", sep="\n")
tmp1 <- gsub(","," ",tmp1)
tmp2 <- strsplit(tmp1, "[[:space:]]+")
names(tmp2) <- sapply(tmp2, `[[`, 1)
annotation <- lapply(tmp2, `[`, -1)
go_anno <- list2df(annotation,col1 = "GO",col2="GENE")
go_anno <- full_join(go_anno, goterms, by = "GO")
go_anno_bp <- go_anno[which(go_anno$ONTOLOGY == "BP"),]
go_anno_mf <- go_anno[which(go_anno$ONTOLOGY == "MF"),]
go_anno_cc <- go_anno[which(go_anno$ONTOLOGY == "CC"),]
go2gene = go_anno[,c("GO","GENE")]
go2name = go_anno[,c("GO","NAMES")]
go2gene_bp = go_anno_bp[,c("GO","GENE")]
go2name_bp = go_anno_bp[,c("GO","NAMES")]
go2gene_mf = go_anno_mf[,c("GO","GENE")]
go2name_mf = go_anno_mf[,c("GO","NAMES")]
go2gene_cc = go_anno_cc[,c("GO","GENE")]
go2name_cc = go_anno_cc[,c("GO","NAMES")]

### Convert GOSeq GO to DOSE GO
filelist = list.files(pattern = "*-UP.subset")
for(file in filelist){
	gene_table <- read.delim(paste0(dir,"/",file), sep="\t",stringsAsFactors=F)
	GOlist_enriched <- paste0(gsub("-UP.subset","-UP.enriched.GO.txt",file))
	GOlist_depleted <- paste0(gsub("-UP.subset","-UP.depleted.GO.txt",file))
	
	if(!file.exists(GOlist_enriched)){ next} else{
		goseq_res <- read.delim(paste0(dir,"/",GOlist_enriched),sep="\t",stringsAsFactors=F)
		cat("\n##############################################################################\n\nWorking with: ", GOlist_enriched,"\n\n##############################################################################\n\n")
		
		goseq_tmp <- goseq_res
		if(dim(goseq_tmp)[1] == 0){next}
		tmp <- go_anno[go_anno$GENE %in% unlist(strsplit(goseq_tmp$gene_ids,", ",fixed=T)),]
		gr_size <- as.character(length(unique(tmp$GENE)))
		bg_size <- as.character(length(which(!is.na(unique(go_anno$GENE)))))
		rownames(goseq_tmp) <- goseq_tmp$category
		goseq_res2 <- goseq_tmp[,c(1,6,4,5,2,8,8,10,4,7)]
		colnames(goseq_res2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue", "geneID","Count","Ontology")
		goseq_res2$GeneRatio <- paste0(goseq_res2$GeneRatio,"/",gr_size)
		goseq_res2$BgRatio <- paste0(goseq_res2$BgRatio,"/",bg_size)
		goseq_res2$geneCode <- goseq_res2$geneID
		goseq_res2$geneID <- gsub(", ","/",goseq_res2$geneID)
		ind<-which(!is.na(unique(go_anno$GENE)))
		universe_genes <- go_anno$GENE[ind]
		
		cat("\nFixing sample genes list...\n")
		gene_table <- rownames(gene_table)
		missing_gtable <- setdiff(gene_table,mapping$V1)
		matching_gtable <- intersect(gene_table,mapping$V1)
		new_gtable <- as.data.frame(matching_gtable)
		colnames(new_gtable) <- "V1"
		new_gtable2 <- left_join(new_gtable, mapping, "V1")
		new_gtable3 <- new_gtable2$V2
		final_gtable <- c(new_gtable3,missing_gtable)
		
		cat("\nFixing gene annotations...\n")
		for(i in 1:length(goseq_res2[,8])){
		gnames1 <- as.data.frame(unlist(strsplit(x=goseq_res2[i,8], split="/")))
		colnames(gnames1) <- "V1"
		gnames2 <- left_join(gnames1, mapping, "V1")
		goseq_res2[i,8] <- paste(unlist(gnames2$V2),collapse="/")
		}
		
		cat("\nFixing universe genes list...\n")
		missing_utable <- setdiff(universe_genes,mapping$V1)
		matching_utable <- intersect(universe_genes,mapping$V1)
		new_utable <- as.data.frame(matching_utable)
		colnames(new_utable) <- "V1"
		new_utable2 <- left_join(new_utable, mapping, "V1")
		new_utable3 <- new_utable2$V2
		final_utable <- c(new_utable3,missing_utable)
		
		geneSets <- as.list(strsplit(goseq_res2$geneID,"/",fixed=T))
		names(geneSets) <- goseq_res2$ID
		
	
		my_conv_GO <- new("enrichResult",
		readable = FALSE,
		result = goseq_res2[which(goseq_res2$Ontology == "BP"),],
		pvalueCutoff = 0.05,
		pAdjustMethod = "BH",
		qvalueCutoff = 0.2,
		organism = "Rhizoprionodon taylori",
		ontology = "BP",
		gene = final_gtable,
		keytype = "UNKNOWN",
		universe = universe_genes,
		gene2Symbol = character(0),
		geneSets = geneSets)

		outfile <- goseq_res2[,c(1,2,10,3,4,5,6,11,8)]
		outfile_bp <- goseq_res2[which(goseq_res2$Ontology == "BP"),]
		outfile_cc <- goseq_res2[which(goseq_res2$Ontology == "CC"),]
		outfile_mf <- goseq_res2[which(goseq_res2$Ontology == "MF"),]
		oldname <- 	paste0(GOlist_enriched)
		outname <- paste0(subDir_enriched_tables,gsub(".enriched.GO.txt",".GO_enrichment_table.txt",basename(oldname)))
		outname_bp <- paste0(subDir_enriched_tables,gsub(".enriched.GO.txt",".GO_enrichment_table_BP.txt",basename(oldname)))
		outname_cc <- paste0(subDir_enriched_tables,gsub(".enriched.GO.txt",".GO_enrichment_table_CC.txt",basename(oldname)))
		outname_mf <- paste0(subDir_enriched_tables,gsub(".enriched.GO.txt",".GO_enrichment_table_MF.txt",basename(oldname)))
		write.table(x=outfile, file=outname, quote=F, sep="\t", row.names=F)
		write.table(x=outfile_bp, file=outname_bp, quote=F, sep="\t", row.names=F)
		write.table(x=outfile_cc, file=outname_cc, quote=F, sep="\t", row.names=F)
		write.table(x=outfile_mf, file=outname_mf, quote=F, sep="\t", row.names=F)
		
		if(!length(grep("BP",goseq_res$ontology)) > 0){next}	
		my_dot <- dotplot(my_conv_GO, showCategory=50)
		plot_name <- paste0(subDir_enriched_plots,gsub(".enriched.GO.txt",".top50_BP_enriched_dotplot.pdf",basename(oldname)))
		ggsave(plot=my_dot, file=paste0(plot_name), device="pdf", width=30,height=25,units="cm")
		my_emap <- emapplot(my_conv_GO, showCategory=50)
		plot_name <- paste0(subDir_enriched_plots,gsub(".enriched.GO.txt",".top50_BP_enriched_emapplot.pdf",basename(oldname)))
		ggsave(plot=my_emap, file=paste0(plot_name), device="pdf", width=30,height=25,units="cm")
		cat("\nYour plots are saved as the following files:\n")
		print(plot_name)
		print(outname)
		}
}
	
for(file in filelist){
	gene_table <- read.delim(paste0(dir,"/",file), sep="\t",stringsAsFactors=F)
	GOlist_enriched <- paste0(gsub("-UP.subset","-UP.enriched.GO.txt",file))
	GOlist_depleted <- paste0(gsub("-UP.subset","-UP.depleted.GO.txt",file))
	
	if(!file.exists(GOlist_depleted)){ next} else{
		goseq_res <- read.delim(paste0(dir,"/",GOlist_depleted),sep="\t",stringsAsFactors=F)
		cat("\n##############################################################################\n\nWorking with: ", GOlist_depleted,"\n\n##############################################################################\n\n")
		
		goseq_tmp <- goseq_res
		if(dim(goseq_tmp)[1] == 0){next}
		tmp <- go_anno[go_anno$GENE %in% unlist(strsplit(goseq_tmp$gene_ids,", ",fixed=T)),]
		gr_size <- as.character(length(unique(tmp$GENE)))
		bg_size <- as.character(length(which(!is.na(unique(go_anno$GENE)))))
		rownames(goseq_tmp) <- goseq_tmp$category
		goseq_res2 <- goseq_tmp[,c(1,6,4,5,3,8,8,10,4,7)]
		colnames(goseq_res2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue", "geneID","Count","Ontology")
		goseq_res2$GeneRatio <- paste0(goseq_res2$GeneRatio,"/",gr_size)
		goseq_res2$BgRatio <- paste0(goseq_res2$BgRatio,"/",bg_size)
		goseq_res2$geneCode <- goseq_res2$geneID
		goseq_res2$geneID <- gsub(", ","/",goseq_res2$geneID)
		ind<-which(!is.na(unique(go_anno$GENE)))
		universe_genes <- go_anno$GENE[ind]
		
		cat("\nFixing sample genes list...\n")
		gene_table <- rownames(gene_table)
		missing_gtable <- setdiff(gene_table,mapping$V1)
		matching_gtable <- intersect(gene_table,mapping$V1)
		new_gtable <- as.data.frame(matching_gtable)
		colnames(new_gtable) <- "V1"
		new_gtable2 <- left_join(new_gtable, mapping, "V1")
		new_gtable3 <- new_gtable2$V2
		final_gtable <- c(new_gtable3,missing_gtable)
		
		cat("\nFixing gene annotations...\n")
		for(i in 1:length(goseq_res2[,8])){
		gnames1 <- as.data.frame(unlist(strsplit(x=goseq_res2[i,8], split="/")))
		colnames(gnames1) <- "V1"
		gnames2 <- left_join(gnames1, mapping, "V1")
		goseq_res2[i,8] <- paste(unlist(gnames2$V2),collapse="/")
		}
		
		cat("\nFixing universe genes list...\n")
		missing_utable <- setdiff(universe_genes,mapping$V1)
		matching_utable <- intersect(universe_genes,mapping$V1)
		new_utable <- as.data.frame(matching_utable)
		colnames(new_utable) <- "V1"
		new_utable2 <- left_join(new_utable, mapping, "V1")
		new_utable3 <- new_utable2$V2
		final_utable <- c(new_utable3,missing_utable)
		
		geneSets <- as.list(strsplit(goseq_res2$geneID,"/",fixed=T))
		names(geneSets) <- goseq_res2$ID
		
	
		my_conv_GO <- new("enrichResult",
		readable = FALSE,
		result = goseq_res2[which(goseq_res2$Ontology == "BP"),],
		pvalueCutoff = 0.05,
		pAdjustMethod = "BH",
		qvalueCutoff = 0.2,
		organism = "Rhizoprionodon taylori",
		ontology = "BP",
		gene = final_gtable,
		keytype = "UNKNOWN",
		universe = universe_genes,
		gene2Symbol = character(0),
		geneSets = geneSets)

		outfile <- goseq_res2[,c(1,2,10,3,4,5,6,11,8)]
		outfile_bp <- goseq_res2[which(goseq_res2$Ontology == "BP"),]
		outfile_cc <- goseq_res2[which(goseq_res2$Ontology == "CC"),]
		outfile_mf <- goseq_res2[which(goseq_res2$Ontology == "MF"),]
		oldname <- 	paste0(GOlist_depleted)
		outname <- paste0(subDir_depleted_tables,gsub(".depleted.GO.txt",".GO_depletion_table.txt",basename(oldname)))
		outname_bp <- paste0(subDir_depleted_tables,gsub(".depleted.GO.txt",".GO_depletion_table_BP.txt",basename(oldname)))
		outname_cc <- paste0(subDir_depleted_tables,gsub(".depleted.GO.txt",".GO_depletion_table_CC.txt",basename(oldname)))
		outname_mf <- paste0(subDir_depleted_tables,gsub(".depleted.GO.txt",".GO_depletion_table_MF.txt",basename(oldname)))
		write.table(x=outfile, file=outname, quote=F, sep="\t", row.names=F)
		write.table(x=outfile_bp, file=outname_bp, quote=F, sep="\t", row.names=F)
		write.table(x=outfile_cc, file=outname_cc, quote=F, sep="\t", row.names=F)
		write.table(x=outfile_mf, file=outname_mf, quote=F, sep="\t", row.names=F)
		
		if(!length(grep("BP",goseq_res$ontology)) > 0){next}	
		my_dot <- dotplot(my_conv_GO, showCategory=50)
		plot_name <- paste0(subDir_depleted_plots,gsub(".depleted.GO.txt",".top50_BP_depleted_dotplot.pdf",basename(oldname)))
		ggsave(plot=my_dot, file=paste0(plot_name), device="pdf", width=30,height=25,units="cm")
		my_emap <- emapplot(my_conv_GO, showCategory=50)
		plot_name <- paste0(subDir_depleted_plots,gsub(".depleted.GO.txt",".top50_BP_depleted_emapplot.pdf",basename(oldname)))
		ggsave(plot=my_emap, file=paste0(plot_name), device="pdf", width=30,height=25,units="cm")
		cat("\nYour plots are saved as the following files:\n")
		print(plot_name)
		print(outname)
		}
}
