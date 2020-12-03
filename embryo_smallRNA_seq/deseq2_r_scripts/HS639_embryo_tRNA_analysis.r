##
##Analysis of project HS639 tRNA
##RNAseq in mice embryos
##Fathers injected with dexamethasone, sperm extracted, IVF performed and embryos collected
##
##Set working directory
  setwd("/home/corcoba/sanger/dseq2/HS639_smallRNA_embryo_sperm3h_sperm2d/tRNAf_analysis/embryo")
  options("width"=200)
  project.name <- "HS639_tRNA_embryo_"
  dir.create(paste(project.name, "r_output", sep=""))
##
##Source necessary functions & libraries
  library("DESeq2")
  library("pheatmap")
  library("ggplot2")
  library("RColorBrewer")
  source("/home/corcoba/sanger/scripts/r/functions/plots/plot_expression_profile_per_sample.rfunction")
  source("/home/corcoba/sanger/scripts/r/functions/plots/cluster_plot_samples.rfunction")
  source("/home/corcoba/sanger/scripts/r/functions/plots/plot_count_matrix.rfunction")
  ##
##1.Build a Deseq2 data object and perform the differential expression analysis
    #read the data from text file (output of featureCounts)
      mydata<-read.table(file="../tRNA_featurecounts.txt", header=TRUE)
      colnames(mydata)
    #Build the matrix of counts rounding the numbers to the nearest integer (requird by deseq functions).
      count.matrix <- round(mydata[,7:15])
      rownames(count.matrix) <- mydata[,1]
      colnames(count.matrix) <- paste("S",c(1:9),sep="")
    # Study design data frames:
      col.data <- data.frame(rep(c("DEX","CON"), times=c(4,5)))
		colnames(col.data)<-c("Treatment")
		rownames(col.data) <- colnames(count.matrix)
    #Vectors of colours to match experimetal design
      col.order <- c(4,2,1,3,5,6,7,8)#my preferred order to extract colours from the palette() function
      legend.col <- palette() [col.order[1:nlevels(col.data[,1])]] #vector of colours for the legend
      sample.colour <- palette() [col.order[as.numeric(col.data[,1])]] #matching vector of colours for each sample
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
      res <- results(p)
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
  # 3.3 Number of genes with 0 reads (per sample)
    colSums(assay(p.data)==0)
  # 3.4 Number of genes with NON 0 reads per sample
    sample.genes <- colSums(assay(p.data)!=0)
  # 3.5 Plot the number of genes and reads per sample
    png(file=paste(project.name , "r_output/", project.name , "plot_reads_and_genes_per_sample.png", sep=""));
	plot(sample.reads,sample.genes, xlab="Number of reads", ylab="Number of genes", type="n", main="No. of reads and genes per sample")
	  text(sample.reads,sample.genes, labels=names(sample.reads), col=sample.colour)
	  legend("topleft", legend=levels(col.data[,1]), col=legend.col, pch=15)
    dev.off()
    #S6 has quite some more reads than any other sample
  # 3.6 Plot the expression patterns of each sample
	  png(file=paste(project.name , "r_output/", project.name , "plot_expression_pattern_per_sample.png", sep=""), width=480, height=960);
	    par(mfrow=c(length(sample.colour),2))
	    par(mar=c(1,2,2,0))
	    par(oma=c(3,3,3,0))
	      for (i in seq(1:ncol(assay(p.data)))) {
		# 3.6.1 Before normalization:
		  barplot(gene.counts[,i], names.arg="", ylim=c(0,2000), border=sample.colour[i], main=colnames(gene.counts)[i])
		# 3.6.2 After normalization:
		  barplot(n.gene.counts[,i], names.arg="", ylim=c(0,2000), border=sample.colour[i], main=colnames(n.gene.counts)[i])
	      }
	      legend("bottomright", legend=levels(col.data[,1]), col=legend.col, pch=15, xpd=TRUE)
	      title(ylab="Counts", xlab="Genes", main="Not Normalized                Normalized", outer=TRUE, cex.main=2, cex.lab=2, line=0.5)
	  dev.off()
##	    
## 5. Visual analysis
  # 5.1 Sample distances
    #A useful first step in an RNA-seq analysis is to assess overall similarity between samples and see how this fits to the expectation from the experiment’s design. 
      #Cluster plot
	png(file=paste(project.name , "r_output/", project.name , "plot_cluster_analysis.png", sep=""));
	  cluster.plot.samples(rld)
	dev.off()
      #There is no clear clustering of the treatment groups.
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
	      legend("bottomright", legend=levels(col.data[,1]), col=legend.col, pch=19)
	  dev.off()  
	 #
	 #Plot all PC combinations
	  png(file=paste(project.name , "r_output/", project.name , "plot_PCAall.png", sep=""))
	    pairs(result$x,col=ifelse(col.data[,1]=="CON", "blue","red"))
	  dev.off()
	  #No clustering visible in any plot
##
##6 Export results
    #If we consider a fraction of 10% false positives acceptable, we can consider all genes with an adjusted p value below  0.1 as significant.
	    resSig <- subset(res, padj < 0.1)
	    write.csv(resSig, file=paste(project.name, "r_output/", project.name, "_pvalues_padj.1.csv", sep=""))
	  #There are a number of genes with low p-values but with the p anjusted set to NA. 
	  #This means that the gene can be analyzed but the mean counts is very low, so it is removed from the analysis by the independent filtering
	  data.frame(res[c(1:34),])
      #Export results for all genes:
	write.csv(res, file=paste(project.name, "r_output/", project.name, "_pvalues_all_genes.csv", sep=""))
      #Export Counts matrices
	write.csv(n.gene.counts, file=paste(project.name, "r_output/", project.name, "_counts_matrix_normalized.csv", sep=""))
	write.csv(gene.counts, file=paste(project.name, "r_output/", project.name, "_counts_matrix_NONnormalized.csv", sep=""))
      
##
##7 Plots
  #7.1 MA-plot
    #Distribution of the model coefs. across all genes
    #The genes with low counts and high variance have been shrunk towards 0 fold change (as explained in the deseq2 paper).
      png(file=paste(project.name , "r_output/", project.name ,"ma_plot_significant_fdr1.png", sep=""));
      #setEPS()
      #postscript(file=paste(project.name, "r_output/", project.name, "ma_plot_fdr1.ps", sep=""));
	plotMA(res, colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l", main=paste(project.name, "DEX FDR=0.1", sep=""))
	#We can label genes easily
	  sig.names <- rownames(resSig)
	  with(res[sig.names, ], {text(baseMean, log2FoldChange, sub("mmu-","",sig.names), pos=4, col="dodgerblue")})
      dev.off()
	
#7.2 Counts Plots
  #Finally we need to visualize the actual count data
  #I will plot all genes that are significant after FDR 0.1 correction and those removed by the independent filtering.
     sign.no <- 34
    #
    #Plot:
      png(file=paste(project.name , "r_output/", project.name ,"normalized_counts_plot_significant_genes.png", sep=""), width=1200, height=600);
	  plot.count.matrix(n.gene.counts, col.data, main="p<0.1 (FDR corrected). Low counts genes shown but removed from analysis by independent filering", stop=sign.no)
      dev.off()
##
##
##Plots for publication
#subset the counts matrix to the significant genes
    sign.counts <- n.gene.counts[rownames(resSig),]
#Create factor with the names of the significant tRNAs removing the last characters to later aggregate the counts by level
    codon <- factor(substr(rownames(sign.counts),1,14))
    table(codon)
#Make data frame aggregating the counts per level (codon)
    sign.counts2 <- data.frame(aggregate(sign.counts, by=list(codon), FUN="sum"), row.names=1)
#Make matrix for pheatmap
    mat2 <- (sign.counts2 - rowMeans(sign.counts2))/apply(sign.counts2,1,sd)
#Import arial fonts to R postscript
    library("extrafont")
#Plot pheatmap
    setEPS()
    postscript(file=paste(project.name , "r_output/", project.name , "plot_pheatmap_fdr1.eps", sep=""), family="Arial")
        pheatmap(mat2, cluster_cols=FALSE, annotation_col=col.data, color = colorRampPalette(rev(brewer.pal(n = 3, name ="RdYlBu")))(100), main="", show_colnames=FALSE)
            library(grid)
            grid.text(label="Z-Score", 0.82,0.98, gp=gpar(cex=0.8))
    dev.off()
#
#Volcano plot
##Volcano plot
    # Create the volcano plot on ggplot
    #remove the genes with NA padj
    rest<-data.frame(res)
    rest$codon <- substr(rownames(rest),6,14)
    rest$thr <- rest$padj > 0.1
    
    #genes to mark in the plot:
    #volcano.genes <- c("Arg-CCT-2", "Gly-GCC-2", "Gly-GCC-3", "Gly-GCC-4", "Gly-CCC-5")
    # label.index<- which(rest$codon %in% volcano.genes)
    #which labels are repeated
    # rest$codon[label.index]
     #remove the indexes of repeated tRNAs
     #label.index <- label.index[-c(2,5,6,7,8,9,10,11)]
     #rest$tlabel <- ""
     #rest$tlabel[label.index] <- rest$codon[label.index]
     #rest<- rest[!is.na(rest$thr),]
    head(rest)
    setEPS()
  postscript(file=paste(project.name , "r_output/", project.name ,"plot_Volcano_fdr1.eps", sep=""), family="Arial", height=4, width=4)
    ggplot(rest, aes(x = rest$log2FoldChange, y = -log10(rest$pvalue), color= rest$thr)) + 
        geom_point() + 
        xlab("log2 fold change") + 
        ylab("-log10 p-value") + 
        theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
              axis.title = element_text(size = rel(1.25))) +
        scale_color_discrete(name = "q-value", labels = c("<0.1", ">0.1")) +
        #geom_text(aes(label=rest$codon),hjust=1,vjust=0, nudge_x=-0.1, show.legend=FALSE, check_overlap=FALSE, size=3) +
        theme_bw()
        dev.off()

