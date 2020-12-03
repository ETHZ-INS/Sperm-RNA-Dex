##
##Analysis of project HS639 miRNA
##RNAseq in mice embryos
##Fathers injected with dexamethasone, sperm extracted, IVF performed and embryos collected
##
##Set working directory
  setwd("/home/corcoba/sanger/dseq2/HS639_smallRNA_embryo_sperm3h_sperm2d/miRNA_analysis/embryo")
  options("width"=200)
  project.name <- "HS639_miRNA_embryo_"
  dir.create(paste(project.name, "r_output", sep=""))
##
##Source necessary functions & libraries
  library("DESeq2")
  library("pheatmap")
  library("RColorBrewer")
  source("/home/corcoba/sanger/scripts/r/functions/plots/plot_expression_profile_per_sample.rfunction")
  source("/home/corcoba/sanger/scripts/r/functions/plots/cluster_plot_samples.rfunction")
  source("/home/corcoba/sanger/scripts/r/functions/plots/plot_count_matrix.rfunction")
##
##1.Build a Deseq2 data object and perform the differential expression analysis
    #read the data from text file (output of featureCounts)
      mydata<-read.table(file="../miRNA_featurecounts.txt", header=TRUE)
      colnames(mydata)
    #Build the matrix of counts rounding the numbers to the nearest integer (requird by deseq functions).
      count.matrix <- round(mydata[,7:15])
      rownames(count.matrix) <- mydata[,1]
      colnames(count.matrix) <- paste(rep(c("DEX","CON"), times=c(4,5)),c(1:9),sep="")
    # Study design data frames:
      col.data <- data.frame(rep(c("DEX","CON"), times=c(4,5)))
		colnames(col.data)<-c("Treatment")
		rownames(col.data) <- colnames(count.matrix)
    #Vector of colours to match experimetal design
      sample.colour <- rep(c("red","blue"), times=c(4,5))
    #Build the deseq data set
      p.data<-DESeqDataSetFromMatrix(countData=count.matrix, colData=col.data, design=~ Treatment)
    #Pre-filtering the data
      #Remove all rows with basically no info (0 or 1 counts)
	nrow(p.data)
	p.data <- p.data[ rowSums(counts(p.data)) > 1, ]
    #Building the results table
    #Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. 
    #If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level.
      p <- DESeq(p.data)
##
## 2.Results
      res <- results(p, alpha = 0.05)
      #order according to p-value
	res <- res[ order(res$pvalue), ]
    #Extract the counts matrices for further analyses, and order the genes according to the p-value attained in the deseq2 analysis
      n.gene.counts <- counts(p, normalized=TRUE)[rownames(res),]
      gene.counts <- counts(p, normalized=FALSE)[rownames(res),]
##
## 3.Transformations
      #For visual representation of the data it is useful to transform them in order to stabilize the variance across the diferent mean counts for each gene
      #2 functions are available, rlog and vst. rlog tends to work best for small samples (n<30)
	rld <- rlog(p.data, blind = TRUE)
	#blind = FALSE means that differences between Tissues and Treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. 
	#The experimental design is not used in the transformation, only in estimating the global amount of variability in the counts. 
	#For a fully unsupervised transformation, set blind=TRUE (is the default).
##
##4.Summary Statistics
  #Study design
    col.data
  # 3.1 Total number of genes entering analysis (with sufficient info in the dataset)
    nrow(p.data)
  # 3.2 Number of RNA transcripts that aligned to the genome (per sample)
    sample.reads <- colSums(assay(p.data))
    sample.reads
  # 3.3 Number of genes with 0 reads (per sample)
    colSums(assay(p.data)==0)
  # 3.4 Number of genes with NON 0 reads per sample
    sample.genes <- colSums(assay(p.data)!=0)
    sample.genes
  # 3.5 Plot the number of genes and reads per sample
    png(file=paste(project.name , "r_output/", project.name , "plot_reads_and_genes_per_sample.png", sep=""));
	plot(sample.reads,sample.genes, xlab="Number of reads", ylab="Number of genes", type="n", main="No. of reads and genes per sample")
	  text(sample.reads,sample.genes, labels=names(sample.reads), col=sample.colour)
	  legend("topleft", legend=c("CON","DEX"), col=c("blue","red"), pch=15)
    dev.off()
  # 3.6 Plot the expression patterns of each sample
	png(file=paste(project.name , "r_output/", project.name , "plot_expression_pattern_per_sample.png", sep=""), width=480, height=960);
	  plot.expression.profile.per.sample(
						gene.counts, 
						n.gene.counts, 
						col.data[,1], 
						ylim=c(0,2000))
	dev.off()
##	    
## 5. Visual analysis
  # 5.1 Sample distances
    #A useful first step in an RNA-seq analysis is to assess overall similarity between samples and see how this fits to the expectation from the experiment’s design. 
    png(file=paste(project.name , "r_output/", project.name , "plot_cluster_analysis.png", sep=""));
	cluster.plot.samples(rld,rld$Treatment)
      dev.off()
      #There is no clear clustering of the treatment groups
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
	      legend("topright", legend=c("CON","DEX"), col=c("blue","red"), pch=19)
	  dev.off()  
	 #
	 #Plot all PC combinations
	  png(file=paste(project.name , "r_output/", project.name , "plot_PCAall.png", sep=""))
	    pairs(result$x,col=ifelse(col.data[,1]=="CON", "blue","red"))
	  dev.off()
	  #No clustering visible in any plot
##
##6 Export results
    #If we consider a fraction of 5% false positives acceptable, we can consider all genes with an adjusted p value below  0.05 as significant.
    #No significant genes were found
    #The following list shows the genes with p-value<0.05 before correction for multiple comparisons
      data.frame(res[c(1:4),])
      
      #Export results for all genes:
	write.csv(res, file=paste(project.name, "r_output/", project.name, "_pvalues_all_genes_padj05.csv", sep=""))
      #Export Counts matrices
	write.csv(n.gene.counts, file=paste(project.name, "r_output/", project.name, "_counts_matrix_normalized.csv", sep=""))
	write.csv(gene.counts, file=paste(project.name, "r_output/", project.name, "_counts_matrix_NONnormalized.csv", sep=""))
      
##
##7 Plots
  #7.1 MA-plot
    #Distribution of the model coefs. across all genes
    #The genes with low counts and high variance have been shrunk towards 0 fold change (as explained in the deseq2 paper).
      png(file=paste(project.name , "r_output/", project.name ,"ma_plot_significant_fdr05.png", sep=""));
      #setEPS()
      #postscript(file=paste(project.name, "r_output/", project.name, "ma_plot_fdr1.ps", sep=""));
	plotMA(res, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main=paste(project.name, "DEX FDR=0.05", sep=""))
	#We can label genes easily
	  #sig.names <- rownames(resSig)
	  #with(res[sig.names, ], {text(baseMean, log2FoldChange, sub("mmu-","",sig.names), pos=4, col="dodgerblue")})
      dev.off()
	
#7.2 Counts Plots
  #Finally we need to visualize the actual count data
  #All gene counts in a single plot
    #
    #Plot all genes
     #Plot Significant genes  (before correction)
     dex05counts <- n.gene.counts[c(1:4),]
      png(file=paste(project.name , "r_output/", project.name ,"normalized_counts_plot_p<0.05_uncorrected.png", sep=""), width=1200, height=300);
	plot.count.matrix(dex05counts,col.data, main="p<0.05 (uncorrected)")
      dev.off()
#7.3 Gene clustering
    mat <- (dex05counts - rowMeans(dex05counts))/apply(dex05counts,1,sd)
png(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_p<0.05_uncorrected.png", sep=""), width=480, height=480)
  pheatmap(mat, cluster_cols=FALSE, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="Heatmap DEX vs CON p<0.05 (uncorrected)")
    library(grid)
    grid.text(label="Z-Score", 0.94,0.98, gp=gpar(cex=0.8))
dev.off()
##
##Plots for publication
#Import arial fonts to R postscript
    library("extrafont")
#MA plot in postscript
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name ,"ma_plot_fdr05.eps", sep=""), family="Arial")
  plotMA(res, ylim = c(-6,6), colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main="", xlab="Mean normalized counts", ylab="Log2 fold change")
dev.off()
