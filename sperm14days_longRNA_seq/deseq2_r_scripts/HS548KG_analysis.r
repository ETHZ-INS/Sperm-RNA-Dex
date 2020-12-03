##
##Analysis of project KG 548
##RNAseq in mice subjected to Dexamethasone injection. Sperm samples
##
##Set working directory
  setwd("/home/corcoba/sanger/dseq2/HS548_longRNA_DEX_FC/KG_samples")
  options("width"=180)
  project.name <- "HS548KG_"
  dir.create(paste(project.name, "r_output", sep=""))
##
##1.Build a Deseq2 data object and perform the differential expression analysis
    #read the data from text file (output of featureCounts)
      mydata<-read.table(file="exon_featurecounts.txt", header=TRUE)
      colnames(mydata)
    #Build the matrix of counts rounding the numbers to the nearest integer (requird by deseq functions).
      count.matrix <- round(mydata[,16:23])
      rownames(count.matrix) <- mydata[,1]
      colnames(count.matrix) <- c("v9","v10","v12","v14","d7","d8","d9","d10")
    # Study design data frames:
      col.data<-data.frame(rep(c("Control","Dexa"),each=4))
		colnames(col.data)<-c("Treatment")
		rownames(col.data) <- colnames(count.matrix)
    #Build the deseq data set
      library("DESeq2")
      p.data<-DESeqDataSetFromMatrix(countData=count.matrix, colData=col.data, design=~ Treatment)
    # Pre-filtering the data
      #Remove all rows with basically no info (0 or 1 counts)
	nrow(p.data)
	  p.data <- p.data[ rowSums(counts(p.data)) > 1, ]
    #Building the results table
    #Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. 
    #If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level.
      p <- DESeq(p.data, fitType="mean")
      res <- results(p, alpha = 0.05)
##
## 2.Transformations
      #For visual representation of the data it is useful to transform them in order to stabilize the variance across the diferent mean counts for each gene
      #2 functions are available, rlog and vst. rlog tends to work best for small samples (n<30)
	rld <- rlog(p.data, blind = TRUE)
	#blind = FALSE means that differences between Tissues and Treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. 
	#The experimental design is not used in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, set blind=TRUE (is the default).
	#We can see the effect of the transform by plotting 2 samples against each other
	par(mfrow=c(1,2))
	  plot(assay(p.data)[,1],assay(p.data)[,2], type="n")
	    text(assay(p.data)[,1],assay(p.data)[,2], labels=seq_along(assay(p.data)[,1]))
	  plot(assay(rld)[,1],assay(rld)[,2],type="n")
	    text(assay(rld)[,1],assay(rld)[,2], labels=seq_along(assay(rld)[,1]))
##
## 3.Summary Statistics
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
	text(sample.reads,sample.genes, labels=names(sample.reads), col=rep(c("blue","red"),each=4))
      dev.off()
    # 3.6 Plot the expression patterns of each sample
	# 3.6.1 Before normalization:
	  png(file=paste(project.name , "r_output/", project.name , "plot_expression_pattern_per_sample_not_normalized.png", sep=""));
	    par(mfcol=c(4,2))
	    par(mar=c(1,2,2,1))
	      for (i in seq(1:ncol(assay(p.data)))) {
		barplot(assay(p.data)[,i], names.arg="", ylim=c(0,2000), border=ifelse(i<=4,"blue","red"), main=paste("Sample ",colnames(p.data[,i]), sep=""))
	      }
	      legend("topright", legend=c("Veh","Dex"),col=c("blue","red"), pch=20)
	  dev.off()
	# 3.6.2 After normalization:
	  png(file=paste(project.name , "r_output/", project.name , "plot_expression_pattern_per_sample_normalized.png", sep=""));
	    par(mfcol=c(4,2))
	    par(mar=c(1,2,2,1))
	      for (i in seq(1:ncol(assay(p.data)))) {
		barplot(counts(p, normalized=TRUE)[,i], names.arg="", ylim=c(0,2000), border=ifelse(i<=4,"blue","red"), main=paste("Sample ",colnames(p.data[,i]), sep=""))
	      }
	      legend("topright", legend=c("Veh","Dex"),col=c("blue","red"), pch=20)
	    dev.off()
	 #Sample v10 seems to have higher couts for many genes.
	 #But after normalization to the sequencing depth, there is not much difference in the overall pattern of expression of the samples.
    #
    #3.7 Types of transcripts present
	#Read the gene biotypes data base from file
	 biotypes<- read.table("../../databases/processed_gene_biotypes.txt", header=TRUE)
	#Select the genes from the database that appear in the sample data
	 sample.data <- data.frame(assay(p.data))
	 biotypes2<- biotypes[biotypes$gene_ID %in% row.names(sample.data),]
	#Eliminate the gene IDs that are not in the new database
	 biotypes2$gene_ID <- factor(biotypes2$gene_ID)
	 length(levels(biotypes2$gene_ID)) #18014 unique genes in the parsed database
	#The database has repeated lines because some areas of the genome generate transcripts with different functions 
	  dim(biotypes2)
	  #How many of the genes in the sample dataset have more than one biotype in the database?
	    occur<-table(biotypes2$gene_ID)
	    reps<-occur[occur>1]
	    dim(reps) #7769 genes have more than one biotype
	    sum(biotypes2$gene_ID %in% row.names(reps)) #Repeated genes take 20135 rows in the parsed database
	    #Which are the biotypes of the genes with repetition?
	      #subset the database to the repeated genes only
		multip <- biotypes2[biotypes2$gene_ID %in% row.names(reps),]
		table(multip$biotype)
	      #Since we are not interested in the some categories, we can regroup them
		biotypes3 <- biotypes2
		levels(biotypes3$biotype) [c(2,4,10,13)] <- "anti-, non-sense, intron or pseudogene"
		levels(biotypes3$biotype) [c(9,10)] <- "protein coding or processed transcript"
	      #And parse the database removing repeated rows
		biotypes3 <- unique(biotypes3)
		dim(biotypes3)
		#How many of the genes in the sample dataset have more than one biotype in the parsed database?
		  occur3<-table(biotypes3$gene_ID)
		  reps3<-occur3[occur3>1]
		  noreps3<-occur3[occur3==1]
		    dim(reps3) #5705 genes have more than one biotype
		    dim(noreps3) #12309 genes have more than one biotype
		#I will use the repeat info to assign biotype to the sample data
		  sample.data$biotype <- factor(rep("multiple biotypes",18014), levels=c(levels(biotypes3$biotype), "multiple biotypes"))
		  sample.data$biotype[which(row.names(sample.data) %in% row.names(noreps3))] <- biotypes3$biotype[which(biotypes3$gene_ID %in% row.names(noreps3))]
		#And create a cross table to summarize the type of transcripts present in each sample
		  cross.table<- matrix(nrow=15, ncol=8, dimnames=list(levels(sample.data$biotype),colnames(sample.data[,1:8])))
		  for (i in 1:8) { cross.table[,i] <- tapply(sample.data[,i],sample.data$biotype, FUN=sum)}
		  cross.table <- cross.table[2:15,]#remove 1st row (NAs)
		#Turn table into %s
		  cross.table.pc <- prop.table(cross.table, margin=2)
		  cross.table.pc <- cross.table.pc*100
		#We can now plot the transcript biotypes per sample
		  png(file=paste(project.name , "r_output/", project.name , "plot_transcript_biotypes_per_sample.png", sep=""), width=600);
		    par(mar=c(5,5,5,18))
		    par(xpd=TRUE)
		      barplot(cross.table.pc, 
			      beside=FALSE, 
			      axisnames=TRUE,
			      axes=FALSE, 
			      horiz=TRUE, 
			      las=2, 
			      col=rainbow(14), 
			      main="Transcript biotypes per sample", 
			      ylab="", 
			      xlab="% of reads", 
			      legend=TRUE, 
			      args.legend = list(x=par("usr")[2]+3, y=par("usr")[4]-0.5, xjust=0))
			axis(1)
		  dev.off()
		#We can also identify the gene with highest read counts per sample
		  apply(assay(p.data), 2, function(x) max(x, na.rm = TRUE))
		  sample.data[which.max(sample.data$v9),]
		  sample.data[which.max(sample.data$v10),]
		  sample.data[which.max(sample.data$v12),]
		  sample.data[which.max(sample.data$v14),]
		  sample.data[which.max(sample.data$d7),]
		  sample.data[which.max(sample.data$d8),]
		  sample.data[which.max(sample.data$d9),]
		  sample.data[which.max(sample.data$d10),]#It is always the same lincRNA
##	    
## 4. Visual analysis
  # 4.1 Sample distances
    #A useful first step in an RNA-seq analysis is to assess overall similarity between samples and see how this fits to the expectation from the experiment’s design. 
      sampleDists <- dist(t(assay(rld))) #calculate the Euclidean distance between samples
    #visualize them in a heatplot.
      library("pheatmap")
      library("RColorBrewer")
      #To plot the distance matrix with the rows/columns arranged by the distances in our distance matrix, we manually provide sampleDists to the clustering_distance argument 
      #of the pheatmap function. Otherwise the pheatmap function would assume that the matrix contains the data values themselves, and would calculate distances between the rows/columns
      #We also specify a color palette using the colorRampPalette function from the RColorBrewer package.
	sampleDistMatrix <- as.matrix( sampleDists )
	rownames(sampleDistMatrix) <- rld$Treatment
	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
      #Cluster plot
	png(file=paste(project.name , "r_output/", project.name , "plot_cluster_analysis.png", sep=""));
	  pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
	dev.off()
      #There is no clear clustering of the treatment groups
  #
  # 4.2 Principal Component Analysis
    #We can build the plot specifying that the color and shape of the points should reflect experimental design.
      pcaData <- plotPCA(rld, intgroup = c( "Treatment"), returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      library("ggplot2")
      #PCA Plot
	png(file=paste(project.name , "r_output/", project.name , "plot_pca_analysis.png", sep=""));
	  ggplot(pcaData, aes(x = PC1, y = PC2)) +
	    geom_point(size =5, aes(color = Treatment)) + 
	    geom_text(label=pcaData$name, size=3) +
	    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
	  coord_fixed()
	dev.off()
      #Sample v10 seems to be taking most of the variance explained by the 1st Principal Component, it may be an outlier
      #Then PC2 seems to separate the 2 groups nicely
##
## 5 Annotating and exporting results
  #The results table contains the Ensembl gene IDs, but alterative names may be more informative
  #Load the annotation packages
    library("AnnotationDbi")
    library("org.Mm.eg.db") #Mus musculus
    columns(org.Mm.eg.db) #to get a list of the available names
  #We can use the mapIds function to add individual columns to our results table.
    res$symbol <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),	#We provide the row names of our results table as a key
                     column="SYMBOL",		#tells the mapIds function which information we want from org.Mm.eg.db
                     keytype="ENSEMBL",		#specify that the keys we give are in the ENSEMBL column of org.Mm.eg.db
                     multiVals="first")		#tells the function what to do if there are multiple possible values for a single input value.
    res$gene.name <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),	
                     column="GENENAME",		
                     keytype="ENSEMBL",		
                     multiVals="first")                 
  #Now the results have the desired external gene IDs:
    resOrdered <- res[order(res$pvalue),]
    head(resOrdered)
  #If we consider a fraction of 5% false positives acceptable, we can consider all genes with an adjusted p value below 5% = 0.05 as significant.
    #We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down- and upregulation.
      resSig <- subset(res, padj < 0.05)
      resSig <- data.frame(resSig[ order(resSig$padj), ])
      write.csv(resSig, file=paste(project.name, "r_output/", project.name, "pvalues_padj.05.csv", sep=""))
    #No genes were differentially expressed between vehicle- and dexamethasone-injected animals
    #The complete list of genes and the respective p values can be found in the HS548KG_all_gene_results.csv file
      write.csv(as.data.frame(resOrdered), file =paste(project.name , "r_output/", project.name ,"all_gene_results.csv" , sep=""))
##
##6 Plots
  #6.1 MA-plot
    #Distribution of the model coefs. across all genes
    #The genes with low counts and high variance have been shrunk towards 0 fold change (as explained in the deseq2 paper).
      png(file=paste(project.name , "r_output/", project.name ,"plot_ma_plot.png", sep=""));
	plotMA(res, ylim = c(-6,6), colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l",  main=paste(project.name, "FDR=0.05", sep=""))
      dev.off()
##
##
	##Until here the DEseq2 analysis as routinely performed.
	##Because of the potential outlier mentioned above, I wanted to look in more detail into the Principal Component Analysis. 
	##The PCA plot shown above is done by a function in the DESeq2 package, but it only uses the 500 genes with highest variance, as a tradeoff for computation speed.
	##I thought it may be worth looking into a PCA of all genes (~19000) and to plot all PC pairs to find potential clusters.
	##It is not very likely that something turns up though (seing the cluster analysis performed above with the sample distances).
##
##7 My own pca_analysis
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
	  png(file=paste(project.name , "r_output/", project.name , "MY_pca_analysis_plot_PC1_PC2.png", sep=""))
	    plot(result$x[,1],result$x[,2], main="Principal Components Analysis", xlab=paste("PC 1 (", percent.var[1], "% Var)",sep=" "), ylab=paste("PC 2 (", percent.var[2], "% Var)",sep=" "), type="n")
	      text(result$x[,1],result$x[,2], labels=rownames(pc.mat), col=ifelse(col.data[,1]=="Control", "blue","red"))
	  dev.off()
	  #The amount of variance necessary to explain Sample v10 is not as big as before now that all genes are considered
	#Plot all PC combinations
	  png(file=paste(project.name , "r_output/", project.name , "MY_pca_analysis_plot_all_PC_pairs.png", sep=""))
	    pairs(result$x, main="All Principal Component Pairs", col=ifelse(col.data[,1]=="Control", "blue","red"))
	  dev.off()
	  	#Even though no clustering is apparent, treated and control groups seem to separate along the PC2.
	    png(file=paste(project.name , "r_output/", project.name , "MY_pca_analysis_plot_PC4_PC2.png", sep=""))
	      plot(result$x[,4],result$x[,2], main="Principal Components Analysis", xlab=paste("PC 4 (", percent.var[4], "% Var)",sep=" "), ylab=paste("PC 2 (", percent.var[2], "% Var)",sep=" "), type="n")
	      text(result$x[,4],result$x[,2], labels=rownames(pc.mat), col=ifelse(col.data[,1]=="Control", "blue","red"))
	    dev.off()
	    #Each PC is constructed as a linear combination of the original variables (genes)
	#The genes with highest influence on the orientation of PC2 (the ones with biggest coefficients in the linear combination) may be the most relevant to separate the treatment groups
	#Nevertheless, with ~18000 genes involved in the analysis, some coefficients may be very high simply by chance.
	#I will plot all PC2 coefficients along theoretical values obtained from a Gaussian distribution with the same mean and sd to detect if there are some abnormally high coefficients in PC2 or they are all just random noise.
		  pc2.loadings <- result$rotation[,2]
		  pc2.loadings <- pc2.loadings[ order(abs(pc2.loadings), decreasing=TRUE) ]
	      #Plot them
		png(file=paste(project.name , "r_output/", project.name , "MY_pca_analysis_plot_PC2_loadings.png", sep=""))
		  plot(seq(1:length(pc2.loadings)),pc2.loadings, main="PC2 Loadings")
		    #To see if the samples deviate from random noise, we can plot a line showing where normally distributed data with same mean and sd would lie
		      gaussian.data <- rnorm(length(pc2.loadings),mean(pc2.loadings),sd(pc2.loadings))
		      points(seq(1:length(pc2.loadings)), gaussian.data[order(abs(gaussian.data), decreasing=TRUE)], pch=20,col="red", cex=0.5)
		      abline(h=max(abs(gaussian.data)), col="green")
		      abline(h=-max(abs(gaussian.data)), col="green")
		      #The green lines mark the most extreme values expected from a Gaussian distribution. 
		      #There may be other smaller values that do not closely follow a Gaussian distribution,
		      #but we are more interested in the extremes because they are coefficients in a linear combination, and the larger ones drive the orientation of the PC6.
		      legend(x="topright", legend=c("PC2 loadings", "Theoretical Gauss distrib."), col=c("black","red"),pch=15)
		 dev.off()
	      #Interestingly, there are some genes with abnormally low (negative) coefficients
	      #We can now select the genes that fall beyond the green lines
		lo<-rownames(data.frame(pc2.loadings[abs(pc2.loadings)>max(abs(gaussian.data))]))
	      #And plot their weight on the PC4 PC2 plot seen above (the length of the line is proportional to the coefficient of the respective gene in PC2
		pcs<-data.frame(cbind(result$x[,4],result$x[,2]))
		loads<-cbind(result$rotation[,4],result$rotation[,2])
		loads<-loads[lo,]
		rownames(loads) <- res$symbol[lo]
		#My own biplot:
		png(file=paste(project.name , "r_output/", project.name , "MY_pca_analysis_plot_PC4_PC2_biplot.png", sep=""))
		  plot(pcs[,1], pcs[,2], main="PC4 PC2 biplot", type="n", xlab=paste("PC 4 (", percent.var[4], "% Var)",sep=" "), ylab=paste("PC 2 (", percent.var[2], "% Var)",sep=" "))
		    text(pcs[,1], pcs[,2], labels=rownames(pcs), col=ifelse(col.data=="Control", "blue", "red"))
		    par(new=TRUE)
		      #One of the axes must be in the same scale than the previous plot, the x in this case.
			#Because the values in the second plot ara much smaller than in the 1st, 
			#we will need a multipliers to scale them up and spread the arrows a bit in the x direction
			multip<-max(abs(pcs[,1]))
			axis3.lim<- round(1/multip, digits=2)#this will help me to scale the 3d axis correctly to the x axis
		      plot(loads[,1], loads[,2], xlim=c(-1,1), type="n",  xlab="", ylab="", axes=FALSE)
			arrows(0,0,multip*loads[,1], loads[,2], col="green", length=0)
			text(multip*loads[,1], loads[,2], labels=rownames(loads), col="green", cex=0.8)
			axis(4, col.axis="green")
			axis(3, col.axis="green",at=seq(-1,1,length.out=length(seq(-axis3.lim,axis3.lim,0.01))), labels=seq(-axis3.lim,axis3.lim,0.01))
		dev.off()
##
##7.1 Counts Plots
  #Finally we need to visualize the actual count data to see if the gene selection has any meaning
  no.rows<-ceiling(length(lo)/4)
  png(file=paste(project.name , "r_output/", project.name , "MY_pca_analysis_plot_PC2_main_genes_counts.png", sep=""), width=960, height=960)
    par(mfrow=c(no.rows,4))
    par(mar=c(2,2,1,1))
    for (i in 1:length(lo)) {
	plotCounts(p, gene = lo[i], intgroup=c("Treatment"), main=ifelse(is.na(res$symbol[lo[i]]), lo[i], res$symbol[lo[i]]))
	}
	dev.off()
##
## 7.2 Exporting my results
##
 #the genes more important for PC2 (see above)
 res2Ordered <- res[lo,]
 res2Ordered <- res2Ordered[order(res2Ordered$pvalue),]
 #The file HS548KG_MY_pca_analysis_PC2_main_genes.csv contains a list of the selected genes from PC2 as described above with their p-values in the deseq2 analysis
  write.csv(as.data.frame(res2Ordered), file =paste(project.name , "r_output/", project.name ,"MY_pca_analysis_PC2_main_genes.csv" , sep="")) 
##
##Gene Ontology analysis
  #To find out if the set of genes selected have a related function, I run a GO enrichment analysis in the web http://geneontology.org/page/go-enrichment-analysis
  #No significant results were reported for either biological process, molecular function or cellular component.
  #But most of the genes (27) are not annotated and were thus not mapped in the GO analysis
##
##
##Plots for publication
#Import arial fonts to R postscript
    library("extrafont")
    #font_import()
    #fonts()
    loadfonts(device = "postscript")
#MA-plot
    setEPS()
    postscript(file=paste(project.name, "r_output/", project.name, "ma_plot_fdr05.eps", sep=""), family="Arial");
	plotMA(res, ylim = c(-6,6), colNonSig = "black", colSig = "red", colLine = "blue", cex=0.8, bty="l",  main="", xlab="Mean normalized counts", ylab="Log2 fold change")
      dev.off()
#Biotypes plot
setEPS()
postscript(file=paste(project.name , "r_output/", project.name , "plot_transcript_biotypes_per_sample.eps", sep=""), family="Arial");
		    par(mar=c(5,5,17,5))
		    par(xpd=TRUE)
		      barplot(cross.table.pc, 
			      beside=FALSE, 
			      axisnames=TRUE,
			      names.arg=paste(rep(c("Con","Dex"), each=4),c(9:12),sep="."), 
			      axes=FALSE, 
			      horiz=TRUE, 
			      las=2, 
			      col=rainbow(15), 
			      main="", 
			      ylab="", 
			      xlab="% of Reads", 
			      legend=TRUE, 
			      width=0.2, 
			      args.legend = list(x=par("usr")[1], y=par("usr")[4], xjust=0, yjust=0))
			axis(1)
		  dev.off()
