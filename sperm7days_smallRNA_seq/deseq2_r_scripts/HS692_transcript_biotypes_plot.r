##Plot of transcript biotypes for project HS692 
##RNAseq in mouse sperm
##	2 x 2 design
##		Treatment: Dexamethasone (DEX) vs Vehicle (CON) injection
##		Time (between injection and sperm collection): 7 days (7D) or 3 hours (3H)
##
##Set working directory
  setwd("/home/corcoba/sanger/dseq2/HS692_DEX_sperm7d_sperm3h/transcript_biotypes")
  options("width"=200)
  project.name <- "HS692_transcript_biotypes_"
  dir.create(paste(project.name, "r_output", sep=""))
##
##read the data from text files (output of featureCounts)
      mydata<-read.table(file="exon_featurecounts.txt", header=TRUE)
      colnames(mydata)
    #Build the matrix of counts rounding the numbers up to the nearest integer (required by deseq functions).
      count.matrix <- round(mydata[,c(13:21,7:12)])
	rownames(count.matrix) <- mydata[,1]
	colnames(count.matrix) <- paste(rep(c("CON","DEX","CON","DEX"), times=c(4,4,4,3)),rep(c("7D","3H"), times=c(8,7)),c(1:15),sep="")
	nrow(count.matrix)
	#remove all rows with nearly no info (0 counts for almost all samples)
	  count.matrix <- count.matrix[ rowSums(count.matrix) > 1, ]
	  head(count.matrix)
## Study design data frames:
      treat <- rep(c("CON","DEX","CON","DEX"), times=c(4,4,4,3))
      time <- rep(c("7D","3H"), times=c(8,7))
      comb <- paste(treat,time, sep="_")
      col.data <- data.frame(cbind(treat,time,comb))
		colnames(col.data)<-c("Treatment","Time","Combined")
		rownames(col.data) <- colnames(count.matrix)
    #Vectors of colours to match experimetal design
      legend.col<- c("blue", "light blue", "red", "orange")
      sample.colour <- legend.col[col.data[,3]]
## Types of transcripts present
  #Read the gene biotypes data base from file
    biotypes<- read.table("/home/corcoba/sanger/dseq2/databases/processed_gene_biotypes.txt", header=TRUE)
  #Select the genes from the database that appear in the sample data
    biotypes2<- biotypes[biotypes$gene_ID %in% row.names(count.matrix),]
  #Eliminate the gene IDs that are not in the new database
	 biotypes2$gene_ID <- factor(biotypes2$gene_ID)
	 paste(length(levels(biotypes2$gene_ID)),"unique genes in the parsed database", sep=" ")
  #The database has repeated lines because some areas of the genome generate transcripts with different functions 
	  dim(biotypes2)
	  #How many of the genes in the sample dataset have more than one biotype in the database?
	    occur<-table(biotypes2$gene_ID)
	    reps<-occur[occur>1]
	    paste(dim(reps), "genes have more than one biotype", sep=" ")
	    sum(biotypes2$gene_ID %in% row.names(reps)) #Repeated genes take 21920 rows in the parsed database
	    #Which are the biotypes of the genes with repetition?
	      #subset the database to the repeated genes only
		multip <- biotypes2[biotypes2$gene_ID %in% row.names(reps),]
		table(multip$biotype)
	      #Since we are not interested in the some categories, we can regroup them
		biotypes3 <- biotypes2
		levels(biotypes3$biotype)
		levels(biotypes3$biotype) [c(2,4,10,13)] <- "anti-, non-sense, intron or pseudogene"
		levels(biotypes3$biotype)
		levels(biotypes3$biotype) [c(9,10)] <- "protein coding or processed transcript"
	      #And parse the database removing repeated rows
		biotypes3 <- unique(biotypes3)
		dim(biotypes3)
		#How many of the genes in the sample dataset have more than one biotype in the parsed database?
		  occur3<-table(biotypes3$gene_ID)
		  reps3<-occur3[occur3>1]
		  noreps3<-occur3[occur3==1]
		     paste(dim(reps3), "genes have more than one biotype", sep=" ")
		     paste(dim(noreps3), "genes have one biotype", sep=" ")
		     multip3 <- biotypes3[biotypes3$gene_ID %in% row.names(reps3),]
		table(multip3$biotype)
##I will use the repeat info to assign biotype to the sample data
    sample.data <- count.matrix
    sample.data$biotype <- factor(rep("multiple biotypes",dim(sample.data)[1]), levels=c(levels(biotypes3$biotype), "multiple biotypes"))
    sample.data$biotype[which(row.names(sample.data) %in% row.names(noreps3))] <- biotypes3$biotype[which(biotypes3$gene_ID %in% row.names(noreps3))]
  #And create a cross table to summarize the type of transcripts present in each sample
    cross.table<- matrix(nrow=15, ncol=15, dimnames=list(levels(sample.data$biotype),colnames(sample.data[,1:15])))
		  for (i in 1:15) { cross.table[,i] <- tapply(sample.data[,i],sample.data$biotype, FUN=sum)}
		  cross.table <- cross.table[c(2:13,15),]#remove 1st & 14th rows (NAs)
		  #Add tRNA data to the cross table
		    tRNAs <- read.csv(file="../tRNA_analysis/HS692_tRNA_r_output/HS692_tRNA_counts_matrix_NONnormalized.csv", row.names=1)
		    cross.table <- rbind(cross.table,colSums(tRNAs))
		    rownames(cross.table)[14] <- "tRNA"
  #Turn table into %s
    cross.table.pc <- prop.table(cross.table, margin=2)
    cross.table.pc <- cross.table.pc*100
##We can now plot the transcript biotypes per sample
		  png(file=paste(project.name , "r_output/", project.name , "plot_transcript_biotypes_per_sample.png", sep=""), width=600);
		    par(mar=c(6,6,5,18))
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
			#rotate 60 degrees, srt=60
			  axis(1)
		  dev.off()