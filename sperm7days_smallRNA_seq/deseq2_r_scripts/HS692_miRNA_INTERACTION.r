##Analysis of project HS692 miRNA INTERACTION
##RNAseq in mouse sperm
##	2 x 2 design
##		Treatment: Dexamethasone (DEX) vs Vehicle (CON) injection
##		Time (between injection and sperm collection): 7 days (7D) or 3 hours (3H)
##
##Set working directory
  setwd("/home/corcoba/sanger/dseq2/HS692_DEX_sperm7d_sperm3h/miRNA_analysis")
  options("width"=200)
  project.name <- "HS692_miRNA_"
  dir.create(paste(project.name, "r_output", sep=""))
##
##Source necessary functions & libraries
  library("DESeq2")
  library("pheatmap")
  library("RColorBrewer")
  source("/home/corcoba/sanger/scripts/r/functions/plots/plot_expression_profile_per_sample.rfunction")
  source("/home/corcoba/sanger/scripts/r/functions/plots/cluster_plot_samples.rfunction")
  ##
##1.Build a Deseq2 data object and perform the differential expression analysis
    #read the data from text files (output of featureCounts)
      mydata<-read.table(file="miRNA_featurecounts.txt", header=TRUE)
      colnames(mydata)
    #Build the matrix of counts rounding the numbers up to the nearest integer (required by deseq functions).
      count.matrix <- round(mydata[,c(13:21,7:12)])
	rownames(count.matrix) <- sub("mmu-","",mydata[,1])
	colnames(count.matrix) <- paste(rep(c("CON","DEX","CON","DEX"), times=c(4,4,4,3)),rep(c("7D","3H"), times=c(8,7)),c(1:15),sep="")
	head(count.matrix)
    # Study design data frames:
      treat <- rep(c("CON","DEX","CON","DEX"), times=c(4,4,4,3))
      time <- rep(c("7D","3H"), times=c(8,7))
      comb <- paste(treat,time, sep="_")
      col.data <- data.frame(cbind(treat,time,comb))
		colnames(col.data)<-c("Treatment","Time","Combined")
		rownames(col.data) <- colnames(count.matrix)
    #Vectors of colours to match experimetal design
      legend.col<- c("blue", "light blue", "red", "orange")
      sample.colour <- legend.col[col.data[,3]]    
    #Build the deseq data set
      p.data<-DESeqDataSetFromMatrix(countData=count.matrix, colData=col.data, design=~ Treatment + Time + Treatment:Time)
    #Pre-filtering the data
      #Remove all rows with basically no info (0 or 1 counts)
	nrow(p.data)
	p.data <- p.data[ rowSums(counts(p.data)) > 1, ]
    #Building the results table
    #Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. 
    #If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level.
      p <- DESeq(p.data)
      resultsNames(p)
##
## 2.Results
      #Interaction:
	# the interaction term, answering: is the condition effect *different* across genotypes? (see ?results "Example 2")
	  res.int <- results(p, name = c("TreatmentDEX.Time7D"), alpha=0.05)
	  res.int <- res.int[ order(res.int$pvalue), ]
      #Treatment:
	  res.treat <- results(p, name= c("Treatment_DEX_vs_CON"), alpha=0.05)
	  res.treat <- res.treat[ order(res.treat$pvalue), ]
      #Time:
	  res.time <- results(p, name= c("Time_7D_vs_3H"), alpha=0.05)
	  res.time <- res.time[ order(res.time$pvalue), ]
      #Extract the counts matrices for further analyses, and order the genes according to the p-value attained in the deseq2 analysis
	  n.gene.counts <- counts(p, normalized=TRUE)
	  gene.counts <- counts(p, normalized=FALSE)
## 3.Transformations
      #For visual representation of the data it is useful to transform them in order to stabilize the variance across the diferent mean counts for each gene
      #2 functions are available, rlog and vst. rlog tends to work best for small samples (n<30)
	rld <- rlog(p.data, blind = TRUE)
##
##4.Summary Statistics
  #Study design
    col.data
  # 4.1 Total number of genes entering analysis (with sufficient info in the dataset)
    nrow(p.data)
  # 4.2 Number of RNA transcripts that aligned to the genome (per sample)
    sample.reads <- colSums(assay(p.data))
	sample.reads
  # 4.3 Number of genes with 0 reads (per sample)
    colSums(assay(p.data)==0)
  # 4.4 Number of genes with NON 0 reads per sample
    sample.genes <- colSums(assay(p.data)!=0)
	sample.genes
  # 4.5 Plot the number of genes and reads per sample
    png(file=paste(project.name , "r_output/", project.name , "plot_reads_and_genes_per_sample.png", sep=""));
	plot(sample.reads,sample.genes, xlab="Number of reads", ylab="Number of genes", type="n", main="No. of reads and genes per sample")
	  text(sample.reads,sample.genes, labels=names(sample.reads), col=sample.colour)
	  legend("topleft", legend=levels(col.data[,3]), col=legend.col, pch=15)
    dev.off()
  # 4.6 Plot the expression patterns of each sample
	  png(file=paste(project.name , "r_output/", project.name , "plot_expression_pattern_per_sample.png", sep=""), width=480, height=960);
	    plot.expression.profile.per.sample(
						gene.counts, 
						n.gene.counts, 
						col.data[,3], 
						ylim=c(0,500),
						legend.colour=legend.col)
	  dev.off()
##	    
## 5. Visual analysis
  # 5.1 Sample distances
    #A useful first step in an RNA-seq analysis is to assess overall similarity between samples and see how this fits to the expectation from the experiment’s design. 
      #Cluster plot
	png(file=paste(project.name , "r_output/", project.name , "plot_cluster_analysis.png", sep=""));
	  cluster.plot.samples(rld, col.data[,3])
	dev.off()
      #No clustering resembling experimental design
      #
  # 5.2 Principal Component Analysis
    	pc.mat<- t(assay(rld))
	result<-prcomp(pc.mat, scale=F)
	#Calculate % of variance explained by PCs
	  ls(result)
	  summary(result)
	  total.var<- sum(result$sdev^2)
	#The total variance explained by the PCs should be equal to the total variance explained by all genes
	  Gene.variances <- apply(pc.mat,2,var) #Are the variances for each gene
	  Total.gene.vars <- sum(Gene.variances) #Should be equal to total.var
	#So the % var explained by each PC is
	  percent.var<-round(result$sdev^2*100/total.var)
	#Plot first 2 PCs
	  png(file=paste(project.name , "r_output/", project.name , "plot_PCA12.png", sep=""))
	    plot(result$x[,1],result$x[,2], main="Principal Component Analysis", xlab=paste("PC 1 (", percent.var[1], "% Var)",sep=" "), ylab=paste("PC 2 (", percent.var[2], "% Var)",sep=" "),type="n", bty="n")
	      rect(xleft=min(result$x-10),ybottom=min(result$x-10), xright=max(result$x+10), ytop=max(result$x+10), col="#eeeeee")
	      grid(col = "white")
	      text(result$x[,1],result$x[,2], labels=rownames(pc.mat), col=sample.colour)
	      legend("bottomright", legend=levels(col.data[,3]), col=legend.col, pch=19)
	  dev.off()  
	 #Plot all PC combinationspng(file=paste(project.name , "r_output/", png(file=paste(project.name , "r_output/", project.name , "plot_expression_pattern_per_sample.png", sep=""), width=480, height=960);
	  png(file=paste(project.name , "r_output/", project.name , "plot_PCAall.png", sep=""), height=720, width=720)
	    pairs(result$x,col=sample.colour)
	  dev.off()  
##6. Export results
    #If we consider a fraction of 5% false positives acceptable, we can consider all genes with an adjusted p value below  0.05 as significant.
    #No significant genes found
    #Genes with p-value<0.05 before correction for MC
      #Interaction:
	data.frame(res.int[c(1:20),])
      #Treatment:
	data.frame(res.treat[c(1:13),])
      #Time:
	data.frame(res.time[c(1:24),])
      #
      #Export results for all genes:
	write.csv(res.int, file=paste(project.name, "r_output/", project.name, "pvalues_INTERACTION_padj05.csv", sep=""))
	write.csv(res.treat, file=paste(project.name, "r_output/", project.name, "pvalues_TREATMENT_padj05.csv", sep=""))
	write.csv(res.time, file=paste(project.name, "r_output/", project.name, "pvalues_TIME_padj05.csv", sep=""))
      #Export Counts matrices
	write.csv(n.gene.counts, file=paste(project.name, "r_output/", project.name, "counts_matrix_normalized.csv", sep=""))
	write.csv(gene.counts, file=paste(project.name, "r_output/", project.name, "counts_matrix_NONnormalized.csv", sep=""))
##7 Plots
  #7.1 MA-plot
    #Distribution of the model coefs. across all genes
    #The genes with low counts and high variance have been shrunk towards 0 fold change (as explained in the deseq2 paper).
    #Interaction:
      png(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_INTERACTION.png", sep=""));
	  plotMA(res.int, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="Interaction FDR=0.05", ylim=c(-8,8))
	  #We can label genes easily
	    #sig.names <- rownames(res[c(1:4),])
	    #with(res[sig.names, ], {text(baseMean, log2FoldChange, sig.names, pos=4, col="dodgerblue")})
      dev.off()
    #Treatment:
      png(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TREATMENT.png", sep=""));
	  plotMA(res.treat, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="Treatment FDR=0.05", ylim=c(-8,8))
      dev.off()
    #Time:
      png(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TIME.png", sep=""));
	  plotMA(res.time, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="Time FDR=0.05", ylim=c(-8,8))
      dev.off()
  #
  #7.2 Counts Plots
  #7.3 Gene clustering  
    anno <- data.frame(col.data[,c(1,2)], row.names=rownames(col.data))
##3. Plot heatmap CON DEX for the genes that were significant in a previous experiment (with injection and sperm collection separated by 14 days)
#List of genes significant in the previous experiment
  sign.14d <- c("miR-30a-5p","miR-30d-5p","miR-30e-5p","miR-872-5p","miR-328-3p","miR-7210-5p","miR-7229-5p","miR-7232-3p","miR-7241-3p")
#Matrix of counts to use:
  icounts <- n.gene.counts[sign.14d,c(1:8,9:15)]
#Pheatmap
  mat <- (icounts - rowMeans(icounts))/apply(icounts,1,sd)
  png(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_comparison_to_14days.png", sep=""), width=480, height=480)
    pheatmap(mat, annotation_col=anno, cluster_cols=FALSE, cluster_rows=FALSE, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="DEX vs CON (sign. at 14 days)")
      library(grid)
      grid.text(label="Z-Score", 0.8,0.98, gp=gpar(cex=0.8))
  dev.off()
##
##Plots for publication
#Import arial fonts to R postscript
    library("extrafont")
    #font_import()
    #fonts()
    loadfonts(device = "postscript")
#MA plots
    #Interaction:
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_INTERACTION.eps", sep=""), family="Arial");
	  plotMA(res.int, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="", ylim=c(-6,6), xlab="Mean of normalized counts", ylab="Log2 fold change")
      dev.off()
    #Treatment:
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TREATMENT.eps", sep=""), family="Arial");
	  plotMA(res.treat, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="", ylim=c(-6,6), xlab="Mean of normalized counts", ylab="Log2 fold change")
      dev.off()
    #Time:
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TIME.eps", sep=""), family="Arial");
	  plotMA(res.time, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="", ylim=c(-6,6), xlab="Mean of normalized counts", ylab="Log2 fold change")
      dev.off()
#Pheatmap of the genes that were significant in the 14 days experiment (this is not the same as above, the 14d experiment was repeated with higher sequencing depth more recently
#obtain the data from the 14 days experiment:
    data14 <- read.csv(file="/home/corcoba/sanger/dseq2/FGCZ6510_dex_free_and_bound/deduplicated_samples/miRNA/FGCZ6510_Free_miRNA_dedup_r_output/FGCZ6510_Free_miRNA_dedup_pvalues_all_genes_padj05.csv", row.names=1)
#select the names of the significant genes in the 14d experiment
    data14$padj
    #The 1st 22 genes have a pval<0.05
    names14 <- rownames(data14[c(1:22),])
#Matrix of counts to use:
    #are all miRNAs sinificant in the 14d experiment present in the current dataset?
        names14 %in% rownames(n.gene.counts)
    #no
        names14 <- names14[-c(5,12)]
#subset the matrix
    counts14 <- n.gene.counts[names14,]
#Pheatmap
  mat14 <- (counts14 - rowMeans(counts14))/apply(counts14,1,sd)
  setEPS()
    postscript(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_comparison_to_14days.eps", sep=""), family="Arial")
    pheatmap(mat14, annotation_col=anno, cluster_cols=FALSE, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="", show_colnames=FALSE)
      library(grid)
      grid.text(label="Z-Score", 0.82,0.98, gp=gpar(cex=0.8))
  dev.off()
