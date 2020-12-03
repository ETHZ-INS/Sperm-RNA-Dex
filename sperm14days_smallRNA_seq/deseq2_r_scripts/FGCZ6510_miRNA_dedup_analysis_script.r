##
##Set working directory
  setwd("/home/corcoba/sanger/dseq2/FGCZ6510_dex_free_and_bound/deduplicated_samples/miRNA")
  options("width"=200)
  project.name <- "FGCZ6510_Free_miRNA_dedup_"
  dir.create(paste(project.name, "r_output", sep=""))
##
##Source necessary functions & libraries
  library("DESeq2")
  library("pheatmap")
  library("RColorBrewer")
  library("ggplot2")
  source("/home/corcoba/sanger/scripts/r/functions/plots/plot_expression_profile_per_sample.rfunction")
  source("/home/corcoba/sanger/scripts/r/functions/plots/cluster_plot_samples.rfunction")
  source("/home/corcoba/sanger/scripts/r/functions/plots/plot_counts_multiple.rfunction")
##
##1.Build a Deseq2 data object and perform the differential expression analysis
    #read the data from text file (output of featureCounts
      mydata<-read.table(file="../miRNA_featurecounts.txt", header=TRUE)
      colnames(mydata)
      #
      #F samples are free RNA
      #B samples are RNA bound to DNA
      #Samples 3 to 10 are Dexamethasone treated
      #Samples 11 to 18 are Controls
      #
    #Build the matrix of counts rounding the numbers to the nearest integer (requird by deseq functions)
      count.matrix <- round(mydata[,c(24:31,32:38,23)])
      colnames(count.matrix)
      rownames(count.matrix) <- sub("mmu-","",mydata[,1])
      colnames(count.matrix) <- paste(rep(c("CON","DEX"), each=8),c(11:18,3:10),sep=".")
    # Study design (write here the experimental group of each sample according to their order in the matrix):
      col.data<-data.frame(rep(c("CON","DEX"), each=8))
		colnames(col.data)<-c("Treatment")
		rownames(col.data) <- colnames(count.matrix)
    #Vectors of colours to match experimetal design
      legend.col<- c("blue","red")
      sample.colour <- legend.col[col.data[,1]]    
    #Build the deseq data set
      p.data<-DESeqDataSetFromMatrix(countData=count.matrix, colData=col.data, design=~ Treatment)
    #Remove all rows with basically no info (0 or 1 counts)
      nrow(p.data)
      p.data <- p.data[ rowSums(counts(p.data)) > 1, ]
    #Results
      p <- DESeq(p.data)
##2.Transformations
  #For visual representation of the data it is useful to transform them in order to stabilize the variance across the diferent mean counts for each gene
  #2 functions are available, rlog and vst. rlog tends to work best for small samples (n<30)
    rld <- rlog(p.data, blind = TRUE)
    #blind = FALSE means that differences between Tissues and Treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. 
    #The experimental design is not used in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, set blind=TRUE (is the default).
##
##3. Summary Statistics
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
	  legend("topleft", legend=c("CON","DEX"), col=legend.col, pch=19)
    dev.off()
  # 3.6 Plot the expression patterns of each sample
	png(file=paste(project.name , "r_output/", project.name , "plot_expression_pattern_per_sample.png", sep=""), width=480, height=1440);
	  plot.expression.profile.per.sample(
						counts(p, normalized=FALSE), 
						counts(p, normalized=TRUE), 
						col.data[,1], 
						ylim=c(0,2000))
	dev.off()
##
##3. Visual analysis
##3.1 Sample distances
   #A useful first step in an RNA-seq analysis issample.colour to assess overall similarity between samples and see how this fits to the expectation from the experiment’s design. 
  #visualize them in a heatplot.
    #Cluster plot
      png(file=paste(project.name , "r_output/", project.name , "plot_cluster_analysis.png", sep=""));
	cluster.plot.samples(rld,rld$Treatment)
      dev.off()
      #There is some clustering of the treatment groups.
##
##3.2 Principal Component Analysis
##
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
	    plot(result$x[,1],result$x[,2], main="Principal Components Analysis", xlab=paste("PC 1 (", percent.var[1], "% Var)",sep=" "), ylab=paste("PC 2 (", percent.var[2], "% Var)",sep=" "),type="n", bty="n")
	      rect(xleft=min(result$x-10),ybottom=min(result$x-10), xright=max(result$x+10), ytop=max(result$x+10), col="#eeeeee")
	      grid(col = "white")
	      text(result$x[,1],result$x[,2], labels=rownames(pc.mat), col=sample.colour)
	      legend("topleft", legend=c("CON","DEX"), col=legend.col, pch=19)
	  dev.off() 
	  #PC1 separates the 2 groups fairly neatly. It explains 26% of the variance
	#Plot all PC combinations
	  png(file=paste(project.name , "r_output/", project.name , "plot_PCAall.png", sep=""), height=960, width=960)
	    pairs(result$x,col=sample.colour)
	  dev.off()  
##
##4. Results
##4.1 Building the results table
    #Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. 
    #If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level.
      res <- results(p, alpha = 0.05)
      res <- res[ order(res$pvalue), ]
    #If we consider a fraction of 5% false positives acceptable, we can consider all genes with an adjusted p value below 0.05 as significant. (22 significant Genes)
      resSig <- subset(res, padj < 0.05)
      data.frame(resSig)
    #Extract the counts matrices for further analyses, and order the genes according to the p-value attained in the deseq2 analysis
      n.gene.counts <- counts(p, normalized=TRUE)[rownames(res),]
      gene.counts <- counts(p, normalized=FALSE)[rownames(res),]
##4.2 Export results
      #For all genes:
	write.csv(res, file=paste(project.name, "r_output/", project.name, "pvalues_all_genes_padj05.csv", sep=""))
      #Export Counts matrices
	write.csv(n.gene.counts, file=paste(project.name, "r_output/", project.name, "counts_matrix_normalized.csv", sep=""))
	write.csv(gene.counts, file=paste(project.name, "r_output/", project.name, "counts_matrix_NONnormalized.csv", sep=""))
##
##5 Plots
##5.1 MA-plot
    #Distribution of the model coefs. across all genes
    #The genes with low counts and high variance have been shrunk towards 0 fold change (as explained in the deseq2 paper).
    png(file=paste(project.name , "r_output/", project.name ,"plot_MA_fdr05.png", sep=""));
        plotMA(res, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main=paste(project.name, "FDR=0.05", sep=""), ylim=c(-3,3))
        #We can label genes easily
            sig.names <- rownames(resSig)
            with(res[sig.names, ], {text(baseMean, log2FoldChange, sig.names, pos=4, col="dodgerblue")})
    dev.off()
       
##
##5.2 Counts Plots
  #subset the counts matrices to the significant genes
    sign.counts <- n.gene.counts[rownames(resSig),]
  #plot
	png(file=paste(project.name , "r_output/", project.name , "plot_gene_counts.png", sep=""), width=960, height=960)
	  plot.counts.multiple(sign.counts, colData(p.data)[,1], sample.colour, normalized=TRUE)
	dev.off()
##5.3 Gene clustering
mat <- sign.counts
mat <- (mat - rowMeans(mat))/apply(mat,1,sd)
png(file=paste(project.name , "r_output/", project.name , "plot_pheatmap.png", sep=""), width=480, height=480)
  pheatmap(mat, cluster_cols=FALSE, annotation_col=col.data, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="Heatmap DEX vs CON (FDR=0.05)")
    library(grid)
    grid.text(label="Z-Score", 0.8,0.99, gp=gpar(cex=0.8))
dev.off()
##
##6 Plots for publication

#Import arial fonts to R postscript
    library("extrafont")
    font_import()
    fonts()
    loadfonts(device = "postscript")
#heatmap plot in postscript
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name , "plot_pheatmap.eps", sep=""), family="Arial")
            pheatmap(mat, cluster_cols=FALSE, annotation_col=col.data, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="", show_colnames=FALSE)
    library(grid)
    grid.text(label="Z-Score", 0.82,0.98, gp=gpar(cex=0.8))
    dev.off()
##6.2 Volcano plot
    # Create the volcano plot on ggplot
    rest<-data.frame(res)
    rest$thr <- rest$padj > 0.05
    rest$thr[is.na(rest$thr)] <- "TRUE"
    setEPS()
  postscript(file=paste(project.name , "r_output/", project.name ,"plot_Volcano_fdr05.eps", sep=""), family="Arial", height=4, width=4)
    ggplot(rest, aes(x = rest$log2FoldChange, y = -log10(rest$pvalue), color= rest$thr)) + 
        geom_point() + 
        xlab("Log2 fold change") + 
        ylab("-Log10 p-value") + 
        theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
              axis.title = element_text(size = rel(1.25))) +
        scale_color_discrete(name = "q-value", labels = c("<0.05", ">0.05")) +
        geom_text(aes(label=ifelse(rest$pvalue<=0.0005,rownames(rest),'')),hjust=1,vjust=0, nudge_x=-0.1, show.legend=FALSE, check_overlap=TRUE, size=3) +
        theme_bw()
        dev.off()
