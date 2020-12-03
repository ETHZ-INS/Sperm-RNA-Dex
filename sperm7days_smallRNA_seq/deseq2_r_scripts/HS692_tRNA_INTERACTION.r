##Analysis of project HS692 tRNA INTERACTION
##RNAseq in mouse sperm
##	2 x 2 design
##		Treatment: Dexamethasone (DEX) vs Vehicle (CON) injection
##		Time (between injection and sperm collection): 7 days (7D) or 3 hours (3H)
##
##Set working directory
  setwd("/home/corcoba/sanger/dseq2/HS692_DEX_sperm7d_sperm3h/tRNA_analysis")
  options("width"=200)
  project.name <- "HS692_tRNA_"
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
      mydata<-read.table(file="tRNA_featurecounts.txt", header=TRUE)
      colnames(mydata)
    #Build the matrix of counts rounding the numbers up to the nearest integer (required by deseq functions).
      count.matrix <- round(mydata[,c(13:21,7:12)])
	rownames(count.matrix) <- mydata[,1]
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
						ylim=c(0,2500),
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
      #Clustering separates very well 7days from 3 hours
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
    #Genes with p-value<0.05
      #Interaction:
	data.frame(res.int[res.int$padj<0.05,])
      #Treatment:
	data.frame(res.int[res.treat$padj<0.05,])
      #Time:
	data.frame(res.int[res.time$padj<0.05,])
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
	  plotMA(res.int, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="Interaction FDR=0.05", ylim=c(-6,6))
	  #We can label genes easily
	    #sig.names <- rownames(res.int[c(1:26),])
	    #with(res.int[sig.names, ], {text(baseMean, log2FoldChange, sig.names, pos=4, col="dodgerblue")})
      dev.off()
    #Treatment:
      png(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TREATMENT.png", sep=""));
	  plotMA(res.treat, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="Treatment FDR=0.05", ylim=c(-6,6))
	  #We can label genes easily
	    #sig.names <- rownames(res.treat[c(1:23),])
	    #with(res.treat[sig.names, ], {text(baseMean, log2FoldChange, sig.names, pos=4, col="dodgerblue")})
      dev.off()
    #Time:
      png(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TIME.png", sep=""));
	  plotMA(res.time, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="Time FDR=0.05", ylim=c(-6,6))
	  #We can label genes easily
	    #sig.names <- rownames(res.time[c(1:190),])
	    #with(res.time[sig.names, ], {text(baseMean, log2FoldChange, sig.names, pos=4, col="dodgerblue")})
      dev.off()
  #
  #7.2 Counts Plots
  #7.3 Gene clustering    
    #List of significant genes for each test
      int.sign <- rownames(res.int[res.int$padj<0.05,])
      treat.sign <- rownames(res.treat[res.treat$padj<0.05,])
      time.sign <- rownames(res.time[res.time$padj<0.05,])
    #Matrices of counts to use:
      int.counts <- n.gene.counts[int.sign,]
      treat.counts <- n.gene.counts[treat.sign,]
      time.counts <- n.gene.counts[time.sign,]
    #Pheatmaps
      int.mat <- (int.counts - rowMeans(int.counts))/apply(int.counts,1,sd)
      treat.mat <- (treat.counts - rowMeans(treat.counts))/apply(treat.counts,1,sd)
      time.mat <- (time.counts - rowMeans(time.counts))/apply(time.counts,1,sd)
      anno <- data.frame(col.data[,c(1,2)], row.names=rownames(col.data))
    #Plots
      png(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_fdr05_INTERACTION.png", sep=""), width=480, height=480)
	pheatmap(int.mat, annotation_col=anno, cluster_cols=FALSE, cluster_rows=FALSE, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="Interaction FDR<0.05")
	library(grid)
	grid.text(label="Z-Score", 0.80,0.98, gp=gpar(cex=0.8))
      dev.off()
      #
      png(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_fdr05_TREATMENT.png", sep=""), width=480, height=480)
	pheatmap(treat.mat, annotation_col=anno, cluster_cols=FALSE, cluster_rows=FALSE, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="Treatment: DEX vs CON FDR<0.05")
	library(grid)
	grid.text(label="Z-Score", 0.80,0.98, gp=gpar(cex=0.8))
      dev.off()
      #
      png(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_fdr05_TIME.png", sep=""), width=480, height=3380)
	pheatmap(time.mat, annotation_col=anno, cluster_cols=FALSE, cluster_rows=FALSE, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="Time: 7 days vs 3 hours FDR<0.05")
	library(grid)
	grid.text(label="Z-Score", 0.80,0.998, gp=gpar(cex=0.8))
      dev.off()
   #7.4 Plot heatmap CON DEX for the genes that were significant in a previous experiment (with injection and sperm collection separated by 14 days)
    #List of genes significant in the previous experiment
      sign.14d <- read.csv(file="/home/corcoba/sanger/dseq2/HS692_DEX_sperm7d_sperm3h/tRNA_analysis/HS692_tRNA_r_output/list_of_significant_tRNASin_14days_experiment.csv")
      sign.14d <- as.vector(sign.14d$x)
    #Matrix of counts to use:
      icounts <- n.gene.counts[sign.14d,c(1:8,9:15)]
    #Pheatmap
      mat <- (icounts - rowMeans(icounts))/apply(icounts,1,sd)
      png(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_comparison_to_14days.png", sep=""), width=480, height=960)
	pheatmap(mat, annotation_col=anno, cluster_cols=FALSE, cluster_rows=FALSE, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="DEX vs CON (sign. at 14 days)")
	library(grid)
	grid.text(label="Z-Score", 0.94,0.98, gp=gpar(cex=0.8))
      dev.off()
##
##Plots for publication
#Import arial fonts to R postscript
    library("extrafont")
    #font_import()
    #fonts()
    #loadfonts(device = "postscript")
#MA plots
    #Interaction:
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_INTERACTION.eps", sep=""), family="Arial");
	  plotMA(res.int, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="", ylim=c(-6,6), xlab="Mean of normalized counts", ylab="Log2 fold change")
	  #We can label genes easily
	    with(res.int["tRNA-Arg-CCT-2-1", ], {text(baseMean, log2FoldChange, "tRNA-Arg-CCT-2-1", pos=4, col="dodgerblue")})
	    with(res.int["tRNA-Arg-CCT-2-1", ], {points(baseMean, log2FoldChange, pch=1, cex=1.5, col="dodgerblue")})
      dev.off()
    #Treatment:
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TREATMENT.eps", sep=""), family="Arial");
	  plotMA(res.treat, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="", ylim=c(-6,6), xlab="Mean of normalized counts", ylab="Log2 fold change")
	  #We can label genes easily
	    with(res.treat["tRNA-Arg-CCT-2-1", ], {text(baseMean, log2FoldChange, "tRNA-Arg-CCT-2-1", pos=4, col="dodgerblue")})
	    with(res.treat["tRNA-Arg-CCT-2-1", ], {points(baseMean, log2FoldChange, pch=1, cex=1.5, col="dodgerblue")})
      dev.off()
    #Time:
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name ,"plot_MAplot_fdr05_TIME.eps", sep=""), family="Arial");
	  plotMA(res.time, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="", ylim=c(-6,6), xlab="Mean of normalized counts", ylab="Log2 fold change")
	  #We can label genes easily
	    with(res.time["tRNA-Arg-CCT-2-1", ], {text(baseMean, log2FoldChange, "tRNA-Arg-CCT-2-1", pos=4, col="dodgerblue")})
	    with(res.time["tRNA-Arg-CCT-2-1", ], {points(baseMean, log2FoldChange, pch=1, cex=1.5, col="dodgerblue")})
      dev.off()
