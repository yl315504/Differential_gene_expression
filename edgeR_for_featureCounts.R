# Rscript edgeR.R $COUNTS_FILE $SAMPLE_FILE $OUT_DIR $EDGER_BCV $CONTRAST_FILE
rm(list=ls())
options(echo=TRUE)
# usage:
usage = "Usage: Rscript script.R  <COUNTS_FILE> <SAMPLE_FILE> <OUT_DIR> <EDGER_BCV> <CONTRAST_FILE> <TPM_file> <AnnotationFile|none>\n"

args_edgeR = commandArgs(trailingOnly = TRUE)
if(length(args_edgeR)!=7) {
  stop("\n    Wrong parameters. \n    ", usage)
}

suppressPackageStartupMessages(library(graphics))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(NMF))


counts_file = args_edgeR[1]
groups_file = args_edgeR[2]
root_output_dir = args_edgeR[3]
edgeR_bcv = as.numeric(args_edgeR[4])
contrast_file = args_edgeR[5]
tpm_file = args_edgeR[6]
anno_file = args_edgeR[7]


counts_file
groups_file
root_output_dir
edgeR_bcv
output_file_prefix = basename(file_path_sans_ext(counts_file))
output_file_prefix
tpm_file
anno_file

topn = 100 # number of top differentially expressed genes to plot in heatmap
image_size_width <- 1024
image_size_height <- 1024

# create output dir
dir.create(file.path(root_output_dir), showWarnings = TRUE)
dir.create(file.path(root_output_dir, "edgeR"), showWarnings = TRUE)
dir.create(file.path(root_output_dir, "edgeR", output_file_prefix), showWarnings = TRUE)
output_dir = paste(root_output_dir,"edgeR",output_file_prefix, sep = "/", collapse = NULL)

output_dir

##########################################
# process counts
counts <- read.delim(counts_file,row.names=1, check.names=FALSE)
keep <- rowSums(cpm(counts)>0) >= 1  # keep only genes with non-zero cpm in at least one sample. in other words, remove genes with all-zero cpms.

#This was for Karam Only
#keep <- rowMeans(counts) >= 10  # keep only genes with non-zero cpm in at least one sample. in other words, remove genes with all-zero cpms.

message("\nFound ", nrow(counts) - sum(keep), " genes with all-zero counts. Keeping the remaining ", sum(keep), ".")
counts <- counts[keep, ]
head(counts)

##########################################
# process sample names, grouping, and contrast information
sample_info = read.delim(groups_file,sep=",", header=TRUE, stringsAsFactors=FALSE)
sample_info = sample_info[, c(1,2)]
rownames(sample_info) = sample_info[, 1]
colnames(sample_info) = c("sampleID", "groupID")
sample_info$groupID = factor(sample_info$groupID)  # convert to factor
design = model.matrix(~0+sample_info$groupID)
colnames(design) = levels(sample_info$groupID)
rownames(design) = sample_info$sampleID
design

counts = counts[ ,rownames(design)]
head(counts)

sink(paste0(output_dir,"/",output_file_prefix,".design.matrix.txt"))
design
sink()

# TODO: check if samples in counts and sample_info are the same.

# read cuffdiff contrast file.
contrast_df <- read.delim(contrast_file, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)  # important to set stringsAsFactors=FALSE
if (nrow(contrast_df) == 0) {
	stop("No comparison specified in ", contrast_file)
}

contrast_list = paste(contrast_df$condition_A, contrast_df$condition_B, sep="-VS-")
contrast_mat <- matrix(0, length(contrast_list), length(levels(sample_info$groupID)))
colnames(contrast_mat) = levels(sample_info$groupID)  # *very* important to set the column names correctly.
rownames(contrast_mat) = contrast_list

for (i in seq_len(nrow(contrast_df))) {
	a = contrast_df$condition_A[i]
	b = contrast_df$condition_B[i]
	if (! all(c(a, b) %in% sample_info$groupID)) {
		stop(contrast_file, " contains an unexpected group name, ", ifelse(a %in% sample_info$groupID, b, a), ". It's not defined in ", groups_file, ": ")
	}
	contrast_mat[i, contrast_df$condition_A[i]] = -1
	contrast_mat[i, contrast_df$condition_B[i]] = 1
}

# each row in contrast_mat specifies one comparison to be analyzed.
contrast_mat

sink(paste0(output_dir,"/",output_file_prefix,".contrast.matrix.txt"))
contrast_mat
sink()

##########################################
# create a DGEList object
d <- edgeR::DGEList(counts=counts, group = sample_info[colnames(counts), "groupID"], lib.size = colSums(counts))

# normalize for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes
# The default method for computing these scale factors uses a trimmed mean of M- values (TMM) between each pair of samples

d <- calcNormFactors(d)

if(grepl(".feature_count.tab",tpm_file)) {
    gene_pos=read.delim(tpm_file,comment.char='#',sep="\t",row.names="Geneid",check.names=FALSE)
    length_pos = gene_pos[c(rownames(counts)),]
    length_stat= cbind(rownames(counts),length_pos$Length)
    write.table(length_stat,file=paste0(output_dir,"/",output_file_prefix,".length.txt"),sep="\t",row.names = F)
    length_file = gene_pos[c(rownames(d$counts)),]
   
    nc = cpm(d)
    col = paste0(colnames(nc), ".CPM")
    row = rownames(nc)
    nc_c = cbind(row,data.frame(nc))
    colnames(nc_c) = c("FEATURE_NAME",col)
    tmp <- paste(output_dir,"/",output_file_prefix,'_', "normalized.counts.CPM.tab",sep="") # output file name for CPM table
    write.table(nc_c, file=tmp, sep="\t",row.names=F)
    
    f=merge(counts,nc_c[,-1],by="row.names")
    colnames(f)[1] <- "FEATURE_NAME"
    
    calc_tpm_from_cpm <- function(x, gene.length) {
        x <- data.matrix(x)
        rpkm.gee <- x*1000/gene.length
        scaling_factor <- colSums(x*1000/gene.length)
        req_tpm= (1000000*rpkm.gee)/scaling_factor[col(rpkm.gee)]
        return (req_tpm)
    }

    tpm <- calc_tpm_from_cpm(nc, gene.length = length_pos$Length)

    tpm_col=paste0(colnames(tpm), ".TPM")
    tpm_row=rownames(tpm)
    tpm = cbind(tpm_row,data.frame(tpm))
    colnames(tpm) = c("FEATURE_NAME",tpm_col)
    tmp_tpm <- paste(output_dir,"/",output_file_prefix,'_', "normalized.counts.TPM.tab",sep="") # output file name for TPM table
    write.table(tpm,file=tmp_tpm,sep="\t",row.names=F)
    
    f=merge(f,tpm,by="FEATURE_NAME")
} else {
    nc = read.delim(tpm_file,row.names=1, check.names=FALSE)
    tpm=nc
    f=merge(counts,nc,by="row.names")
    colnames(f)[1] <- "FEATURE_NAME"
    col = colnames(nc)
    row = rownames(nc)
    nc_c = cbind(row,data.frame(nc))
    colnames(nc_c) = c("FEATURE_NAME",col)
    tmp <- paste(output_dir,"/",output_file_prefix,'_', "normalized.counts.TPM.tab",sep="") # output file name for TPM table
    write.table(nc_c, file=tmp, sep="\t",row.names=F, quote=F)
}

## Draw dendogram for normalized counts
dataset1 <- t(nc_c[,-1])  # ommit FEATURE_NAME
distance_matrix <- dist(dataset1, method = "euclidean") # distance matrix
fit_dm <- hclust(distance_matrix, method="ward")
png(filename = paste(output_dir,"/",output_file_prefix,".normalized.dendogram.png",sep=""), width = image_size_width, height = image_size_height)
plot(fit_dm) # display dendogram
dev.off()

## clustering using  Multidimensional scaling methods
mycols = rainbow(length(levels(d$samples$group )))
if (length(levels(d$samples$group )) <= 8 ) {
	mycols = brewer.pal(n=length(levels(d$samples$group)), name = "Dark2")
}
png(filename = paste(output_dir,"/",output_file_prefix,".sample.clustering_top10percentgenes.png",sep=""), width = image_size_width, height = image_size_height)
plotMDS(d,top=dim(d$counts)[1]/100*10, gene.selection="pairwise",method="bcv", col=mycols[d$samples$group],cex=2 )
dev.off()

png(filename = paste(output_dir,"/",output_file_prefix,".sample.clustering_allgenes.png",sep=""), width = image_size_width, height = image_size_height)
plotMDS(d,top=dim(d$counts)[1], gene.selection = "pairwise",method="bcv", col=mycols[d$samples$group],cex=2 )
dev.off()

##########################################
# differential expression analysis

# Fit a negative binomial generalized log-linear model to the read counts for each gene.
if (edgeR_bcv == 2){ # condition where there are replicates
	message("\nBIOLOGICAL Replicates mode!")

	# estimate common dispersion
	# use the quantile-adjusted conditional maximum likelihood (qCML) method for experiments with single factor
	# The easiest way to share information between genes is to assume that all genes have the same mean-variance relationship, in other words, the dispersion is the same for all the genes
	# However, the truth is that the gene expression levels have non-identical and dependent distribution between genes, which makes the above assumptions too naive
	# A more general approach that allows genewise variance functions with empirical Bayes shrinkage was introduced several years ago (Robinson and Smyth, 2007)
	# Only when using tagwise dispersion will genes that are consistent between replicates be ranked more highly than genes that are not
	# Therefore, the tagwise dispersions are strongly recommended in model fitting and testing for differential expression

	d <- estimateGLMCommonDisp(d,design)
	d <- estimateGLMTrendedDisp(d,design)
	d <- estimateGLMTagwiseDisp(d,design)

	fit <- glmFit(d,design)
} else {  # condition where there are no replicates
  print ("NO Replicates or technical replicates mode!")
  
  if (edgeR_bcv == 0.01){
    print ("RNASeq replicates: Technical")
  }else if (edgeR_bcv == 0.1){
    print ("RNASeq replicates: None (inbreeding organism i.e. mouse, fly etc..)")
  }else if (edgeR_bcv == 0.4){
    print ("RNASeq replicates: None (outbreeding organism i.e. human)")
  }
  
  bcv <- edgeR_bcv
  fit <- glmFit(d,design,bcv^2)
}

# test for differentially expressed genes
# columns for the log2-fold-change, logFC, the average log2-counts-per-million, logCPM, and the two-sided p-value PValue

for (i in seq_along(contrast_list)) {
  message("\nRunning comparision: ", contrast_list[i])
  lrt = glmLRT(fit, contrast=contrast_mat[contrast_list[i],])
  tab <- topTags(lrt, n = Inf)
  my_tab = tab
  tmp <- paste(output_dir,"/",output_file_prefix,'_', contrast_list[i], ".glmLRT.txt",sep="")
  
  col = paste0(contrast_list[i],".",colnames(my_tab))
  row = rownames(my_tab)
  my_tab = cbind(row,data.frame(my_tab))
  colnames(my_tab) = c("FEATURE_NAME",col)
  write.table(my_tab, file=tmp, sep="\t",row.names=F, quote=F)
  f=merge(f,my_tab,by="FEATURE_NAME")
  
  col = colnames(tab)
  row = rownames(tab)
  tab = cbind(row,data.frame(tab))
  colnames(tab) = c("FEATURE_NAME",col)

  # generate a plot of the tagwise log-fold-changes against log-cpm (analogous to an MA-plot for microarray data). DE tags are highlighted on the plot
  # The horizontal blue lines show 4-fold changes
  message("\nGenerating MA plot")
  detags = rownames(tab)[tab$FDR <= 0.05]  # later version of edgeR may use detags <- rownames(topTags(lrt, p.value=0.05))
  tmp <- paste(output_dir,"/", output_file_prefix, '_', contrast_list[i],".MAplot.png",sep="")
  png(tmp, width = image_size_width, height = image_size_height)
  plotSmear(lrt, de.tags=detags, main=contrast_list[i], pch=19, cex=0.5)
  abline(h = c(-2, 2), col = "blue")
  dev.off()
  
  # volcano plot
  message("\nGenerating volcano plot")
  tmp <- paste(output_dir,"/", output_file_prefix, '_', contrast_list[i],".volcano.png",sep="")
  png(tmp, width = image_size_width, height = image_size_height)
  plot(tab$logFC, -log10(tab$PValue), pch=20, cex=1, main=contrast_list[i],
       ylab="-log10PValue", xlab="log2FC", col=adjustcolor(as.numeric(rownames(tab) %in% detags)+1, alpha.f = 0.3))
  abline(v=c(-2, 2), lty=2)
  mtext(text="genes with FDR < 0.05 are highlighted", side=3)
  dev.off()

	#heatmap
	message("\nGenerating heatmap with ", topn, " genes.")
	mat = log2(tpm[rownames(tab)[1:topn], paste0(rownames(d$samples),".TPM")[order(d$samples$group)]] + 0.125)
	head(mat)
	if (any(is.na(cor(t(mat))))) {
		message("Skipping heatmap because some genes have zero standard deviation")
	} else {
		nmf.options(grid.patch=TRUE)
		tmp <- paste(output_dir,"/", output_file_prefix, '_', contrast_list[i],".heatmap.pdf",sep="")
		pdf(tmp,width=10)
		mat.heat=data.matrix(mat)
		par(cex.main=.8)
		heatmap.2(mat.heat,main=paste(contrast_list[i], "Top 100 diff. expr. genes"),
				trace="none",keysize=1,
				Colv=FALSE,dendrogram="row",scale = "row",
				distfun=function(x) as.dist(1 - cor(t(x), method='pearson')),
				col=bluered(75),labRow=rownames(mat),margins=c(10,20),cexCol = 0.7,cexRow = .1
			)
		dev.off()
	}
}

if(anno_file != "none")
{
	anno <- read.delim(anno_file,row.names=1, check.names=FALSE)
	f=merge(f,anno,by="FEATURE_NAME")
}

tmp <- paste(output_dir,"/",output_file_prefix,'.', "combined.txt",sep="") # output all combined results.
write.table(f, row.names=FALSE, file = tmp, sep = "\t", quote = F)