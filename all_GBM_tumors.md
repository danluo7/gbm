# gbm all GBM tumors (049, 011, 024,UNC)

## this continues from raw data processing. the following is Ballgown analysis for differential expression and beyond. 

combining all tumors (replicates merged)


	mkdir -p /home/daniel/ubuntu/workspace/gbm_merged/all_tumors_merged/de/ballgown/ref_only/
	cd /home/daniel/ubuntu/workspace/gbm_merged/all_tumors_merged/de/ballgown/ref_only/

Generate a header file to load into R. The order is the order that the samples get loaded into R object/table. This should also match the shortnames later. when ordering these, try to group similar columns together ie all the tissues togther, so it's easy to compare via c(1:3) to select 3 columns in R. 

this is for ALL samples for principal component analysis (the simplest form of multidimentional scaling), and also a file for pairwise comparisons. since we have a ton of comparisisons, might just not do this for now and only do the PCA. 

file for all samples for PCA: (this is how the script should look like (without the enters inbetween each line):


printf "\"ids\",\"type\",\"path
\"\n\"011_s\",\"011_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_s
\"\n\"011_o\",\"011_organoid\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_o
\"\n\"011_t\",\"011_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_t
\"\n\"011_v\",\"011_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_v

\"\n\"024_s\",\"024_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_s
\"\n\"024_o\",\"024_organoid\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_o
\"\n\"024_t\",\"024_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_t
\"\n\"024_v\",\"024_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_v


\"\n\"049_n\",\"049_neurosphere\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_n
\"\n\"049_s\",\"049_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_s
\"\n\"049_o\",\"049_organoid\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_o
\"\n\"049_t\",\"049_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_t
\"\n\"049_v\",\"049_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_v


\"\n\"unc_s\",\"unc_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_s
\"\n\"unc_t\",\"unc_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_t
\"\n\"unc_v\",\"unc_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_v

\"\n" > all.csv




printf "\"ids\",\"type\",\"path\"\n\"011_s\",\"011_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_s\"\n\"011_o\",\"011_organoid\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_o\"\n\"011_t\",\"011_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_t\"\n\"011_v\",\"011_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_v\"\n\"024_s\",\"024_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_s\"\n\"024_o\",\"024_organoid\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_o\"\n\"024_t\",\"024_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_t\"\n\"024_v\",\"024_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_v\"\n\"049_n\",\"049_neurosphere\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_n\"\n\"049_s\",\"049_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_s\"\n\"049_o\",\"049_organoid\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_o\"\n\"049_t\",\"049_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_t\"\n\"049_v\",\"049_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_v\"\n\"unc_s\",\"unc_slice\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_s\"\n\"unc_t\",\"unc_tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_t\"\n\"unc_v\",\"unc_invitro\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_v\"\n" > all.csv






	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)
	
	dir()


	pheno_data = read.csv("all.csv")  


	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg


	bg_table = texpr(bg, 'all')
	head(bg_table)

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)


	save(bg, file='bg.rda')

	bg




#Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
	
	load('bg.rda')
	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])



	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)



	colnames(gene_expression)


	
	#####row.names(gene_expression)



	dim(gene_expression)


#The data frame is always the format table[rows 1 through whatever, c(columns 1 through whatever,any other column you need added]). To view the first 3 rows of gene #expression data for the first 4 samples plus the sample in the 6th column (1:3 means 1 through 3).
	
	gene_expression[1:3, c(1:4,6)]

#To view gene expression for a single gene by pulling out row names that matches (or ==) BRD4, across all columns (that weird [i,] which means i want every column, saves #typing c(1:10) by just not typing anything)


	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]




	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)
	
	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count



#Many genes will have only 1 transcript, some genes will have several transcripts. Use the 'table()' command #to count the number of times each gene symbol occurs (i.e. #the # of transcripts that have each gene symbol). Then use the 'hist' command to create a #histogram of these counts


## Plot the number of transcripts per gene. 

	pdf(file="all_R_output.pdf") #open up a pdf writer first

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)


## Plot the distribution of transcript sizes as a histogram. lengths will be those of known transcripts. 
#Good QC step: we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts.




	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

	
	min(gene_expression[,"FPKM.1"])
	max(gene_expression[,"FPKM.2"])


	min_nonzero=1


	data_columns=c(1:16)
	short_names=c("011_s","011_o","011_t","011_v","024_s","024_o","024_t","024_v","049_n","049_s","049_t","049_o","049_v","unc_s","unc_t","unc_v")




	#to see all available colors use: colors()
		data_colors=c("tomato1","royalblue1","seagreen1","grey1","tomato2","royalblue2","seagreen2","grey2","black","tomato3","seagreen3","royalblue3","grey3","tomato4","seagreen4","grey4")






	boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 16 samples")



	dev.off() #close off the pdf and start a new one for MDS
	






	
	pdf(file="all_R_output_MDS.pdf")



## Compare the correlation distance between all replicates


#Scree plot: Determine the amount of variance coming from each principal component in a table:

	pc <- princomp(gene_expression[,data_columns],cor=TRUE,scores=TRUE)
	summary(pc)
	plot(pc,type='lines')
	

#Calculate the FPKM sum for all 16 libraries

	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Filter out genes with a grand sum FPKM of less than 10

	i = which(gene_expression[,"sum"] > 10)


#Calculate the correlation between all pairs of data

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r
	



## Plot MDS.
#Convert correlation to distance, and use 'multi-dimensional scaling' to plot the relative differences between libraries, by calculating 2-dimensional #coordinates to #plot points for each library using eigenvectors (eig=TRUE). d, k=2 means 2 dimensions
	
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (PCA) for for 011, 024, 049, unc tumors (all non-zero genes)", xlim=c(-0.4,0.9), ylim=c(-0.4,0.6))
	points(mds$points[,1], mds$points[,2], col=data_colors, cex=1, pch=21)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)


	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (PCA) for for 011, 024, 049, unc tumors (all non-zero genes)", xlim=c(-0.3,0.8), ylim=c(-0.3,0.5))
	points(mds$points[,1], mds$points[,2], col=data_colors, cex=1, pch=20)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)


	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (PCA) for for 011, 024, 049, unc tumors (all non-zero genes)", xlim=c(-0.2,0.2), ylim=c(-0.2,0.2))
	points(mds$points[,1], mds$points[,2], col=data_colors, cex=1, pch=20)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)


#no text labels
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (PCA) for for 011, 024, 049, unc tumors (all non-zero genes)", xlim=c(-0.3,0.8), ylim=c(-0.3,0.5))
	points(mds$points[,1], mds$points[,2], col=data_colors, cex=2, pch=20)
	text(mds$points[,1], mds$points[,2], short_names, col=white)


	dev.off()
	q()
	n
	








## comparing tissue to non-tissue, using stattest function, and output a heatmap


#Need to make a new header file and recode the "type" header,since this will determine what gets compared in a stattest. the stattest function will use type as a covariate and use fpkm as a meansurement. since this function can't compare multiple things, need to make another file called GBM049_all_stattest.csv and make the "type" tissue vs. non-tissue. Then heatmaps can be done comparing invitro to everything else. 

	mkdir -p /home/daniel/ubuntu/workspace/gbm_merged/all_tumors_merged/de/ballgown/ref_only/stattest
	cd /home/daniel/ubuntu/workspace/gbm_merged/all_tumors_merged/de/ballgown/ref_only/stattest




printf "\"ids\",\"type\",\"path
\"\n\"011_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_t
\"\n\"024_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_t
\"\n\"049_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_t
\"\n\"unc_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_t

\"\n\"011_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_s
\"\n\"024_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_s
\"\n\"049_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_s
\"\n\"unc_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_s

\"\n\"011_o\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_o
\"\n\"024_o\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_o
\"\n\"049_o\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_o

\"\n\"011_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_v
\"\n\"024_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_v
\"\n\"049_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_v
\"\n\"unc_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_v

\"\n\"049_n\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_n


\"\n" > all_stattest.csv





printf "\"ids\",\"type\",\"path\"\n\"011_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_t\"\n\"024_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_t\"\n\"049_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_t\"\n\"unc_t\",\"tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_t\"\n\"011_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_s\"\n\"024_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_s\"\n\"049_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_s\"\n\"unc_s\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_s\"\n\"011_o\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_o\"\n\"024_o\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_o\"\n\"049_o\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_o\"\n\"011_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/011_merged/expression/stringtie/ref_only/011_v\"\n\"024_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/024_merged/expression/stringtie/ref_only/024_v\"\n\"049_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_v\"\n\"unc_v\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/unc_merged/expression/stringtie/ref_only/unc_v\"\n\"049_n\",\"non-tissue\",\"/home/daniel/ubuntu/workspace/gbm_merged/049_merged/expression/stringtie/ref_only/049_n\"\n" > all_stattest.csv





#now rerun all the R scripts to see if the resulting MDS looks weird (it might) and see if stattest and heat map work. 


	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("all_stattest.csv")  

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg

	bg_table = texpr(bg, 'all')

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

	save(bg, file='bg.rda')

	bg

	pdf(file="all_R_output_stattest.pdf")


	dir()
		
	load('bg.rda')
		
	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

	colnames(gene_expression)
	dim(gene_expression)

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]

	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)
		
	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")
		
	min(gene_expression[,"FPKM.1"])
	max(gene_expression[,"FPKM.2"])
		
	min_nonzero=1

#make sure the columns AND short names are correct, this is for the names of columns in the heatmap
	data_columns=c(1:16)
short_names=c("011_t","024_t","049_t","unc_t","011_s","024_s","049_s","unc_s","011_o","024_o","049_o","011_v","024_v","049_v","unc_v","049_n")




	#colors()
	
	data_colors=c("green","green","green","green","red","red","red","red","blue","blue","blue","grey","grey","grey","grey","black")



	boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 16 samples")

	



	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
	
	i = which(gene_expression[,"sum"] > 10)

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r




## Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))


## View the distribution of differential expression values as a histogram

#Display only those that are significant according to Ballgown

	sig=which(results_genes$pval<0.05)
	results_genes[,"de"] = log2(results_genes[,"fc"])
	hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) tissue vs non-tissue", main="Distribution of differential expression values")
	abline(v=-2, col="black", lwd=2, lty=2)
	abline(v=2, col="black", lwd=2, lty=2)
	legend("topleft", "Fold-change > 4", lwd=2, lty=2)



#Display the grand expression values from tissue vs invitro and mark those that are significantly differentially expressed. Make sure all the columns are correctly indicated. 

	gene_expression[,"tissue"]=apply(gene_expression[,c(1:4)], 1, mean)
	gene_expression[,"invitro"]=apply(gene_expression[,c(12:15)], 1, mean)

	x=log2(gene_expression[,"tissue"]+min_nonzero)
	y=log2(gene_expression[,"invitro"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="FPKM (log2)", ylab="FPKM (log2)", main="tissue vs invitro FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)



Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
	
	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


Write a simple table of differentially expressed transcripts to an output file

Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]


	o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE) #Order the output by or p-value and then break ties using fold-change

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_R_ballgown_tissue_vs_invitro.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]



#Now for slices vs invitro: same thing, display the grand expression values from slice vs invitro and mark those that are significantly differentially expressed. Make sure all the columns are correctly indicated. 

	gene_expression[,"slice"]=apply(gene_expression[,c(5:8)], 1, mean)
	gene_expression[,"invitro"]=apply(gene_expression[,c(12:15)], 1, mean)

	x=log2(gene_expression[,"slice"]+min_nonzero)
	y=log2(gene_expression[,"invitro"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="FPKM (log2)", ylab="FPKM (log2)", main="slice vs invitro FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)



#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
	
	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


#Write a simple table of differentially expressed transcripts to an output file

#Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]


	o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE) #Order the output by or p-value and then break ties using fold-change

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_R_ballgown_slice_vs_invitro.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]






#### Plot #11 - Create a heatmap to vizualize expression differences between the 16 samples

#Define custom dist and hclust functions for use with heatmaps

	mydist=function(c) {dist(c,method="euclidian")}
	myclust=function(c) {hclust(c,method="average")}

	main_title="sig DE Transcripts for 4 GBM tumors"
	par(cex.main=0.8)
	sig_genes_de=sig_tn_de[,"id"]
	sig_gene_names_de=sig_tn_de[,"gene_name"]

	data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
	heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,labCol=short_names,col=rev(heat.colors(75)))


	dev.off()
	
	
	
