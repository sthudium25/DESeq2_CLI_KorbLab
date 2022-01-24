library("DESeq2")
library("edgeR")
library("tximport")
library("readr")
library(GenomicFeatures)
library(tidyverse)
library(ggrepel)
library(fs)
library(argparse)

today=Sys.Date()
parser = ArgumentParser()
parser$add_argument('-o', '--organism', 
                    default='mouse', 
                    type='str', 
                    help='The organism from which your sequencing files originated [default %(default)s]')
parser$add_argument('-a', '--analysis_name',
                    default=sprintf('Deseq2_%Y%m%d', today),
                    help='A name for the analysis being run [default %(default)s')
parser$add_argument('-p', '--path',
                    help="Full path to your Deseq template file")

args = parser$parse_args()
if (args$organism == 'mouse' ) {
  library(org.Mm.eg.db)
  annots = dir_ls(path = '/Users/Sammy/Desktop/', recurse = TRUE, glob = '*vM19.sqlite', type = 'file')
  annots_filtered = annots[basename(annots) == 'gencode.vM19.sqlite']
  if (length(annots) == 0) {
    download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.annotation.gtf.gz", 
                  "gencode.vM19.annotation.gtf.gz")
    txdb <- makeTxDbFromGFF("gencode.vM19.annotation.gtf.gz")
    saveDb(txdb_ms, file="/Users/Sammy/Desktop/gencode.vM19.sqlite")
  } else {
    txdb_ms <- loadDb(annots[1]) 
  }
} 

if (args$organism == 'human') {
  
  # TODO: repeat for human 
}
  
# This will then create a deseq folder if it isn't already there

folders = dir()
if (!('deseq' %in% folders)) {
  dir.create('deseq')
}
print(getwd())

# Give your current analysis a name. This will be used for the naming of output files
analysis_name = args$analysis_name

columns(txdb)
k <- keys(txdb, "GENEID")
res <- AnnotationDbi::select(txdb, k, "TXNAME", "GENEID") #for every gene, tell me the transcripts that are associated with it
tx2gene <- res[,2:1] #this will show a list that has one column with the transcript name, and another column with the corresponding geneID
head(tx2gene)

sampleInfo = read_csv(file = 'Deseq_CLI_test.csv', col_names = T) # Change to args$path

# This should then point to the quants folder located inside of the current experiment
dir <- file.path("quants/")
list.files(dir)
samplenames <- list.files(dir) #this is the directory with my Salmon outputs and I wirte their names into a vector 
samplenames
samplenames_reorder <- samplenames[c(1:6)] #from the list of files, I select the ones I want for the untreated and LPs comparison (from the list of file names I selected above), and put them in the correct order
samplenames_reorder

# for differential expression between control and one condition alone 
samplenames_subset <- samplenames_reorder[c(1,3,4:6)] #from the list of files, I select the ones I want for the LPS alone and LPS-RMM comparison (from the list of file names I selected above), and put them in the correct order
samplenames_subset
#Write out the treatment conditions: (the levels tells the analysis which conditions to set as a control to compare to)
#Treatment <- factor(c("Control", "Control", "Control", "Dot1Li_4h","Dot1Li_4h","Dot1Li_4h"), levels=c("Control","Dot1Li_4h"))
#Alternatively can 'repeat' a name to avoid having to retype:
Treatment <- factor(c(rep("Ctrl",2),rep("Gtf2iKO",3)), levels=c("Ctrl","Gtf2iKO"))
Treatment
colData <- data.frame(samplenames_subset, Treatment)
colData #this is a table that has the exact file names of the Salmon output, and the other columns are descriptors of the experimental design that will be important for the DESeq2 analysis later on

#Now we can build a vector which points to our quantification files using this column of coldata. We use names to name this vector with the run IDs as well.
files <- file.path(dir,colData$samplenames_subset,"quant.sf")
names(files) <- colData$samplenames_subset
head(files,12)

#use Tximport to read in files correctly. dim gives dimension readout. Should be the number of lines and the number of samples
txi <- tximport(files, type="salmon", tx2gene=tx2gene_ms,  ignoreAfterBar = TRUE)
dim(txi$abundance)
dim(txi$counts)
dim(txi$length)

#Now, we will build a DESeqDataSet from the matrices in tx

dds <- DESeqDataSetFromTximport(txi, colData, ~ Treatment)
dds
#My favorite of these transformation is the vst, mostly because it is very fast, and provides transformed (nearly log-scale) data which is robust to many problems associated with log-transformed data (for more details, see the DESeq2 workflow or vignette ).
#blind=FALSE refers to the fact that we will use the design in estimating the global scale of biological variability, but not directly in the transformation:
vst <- vst(dds, blind=FALSE)
pdf(sprintf('deseq/%s_PCA.pdf', analysis_name))
p1 = plotPCA(vst, "Treatment")
p1 + ggtitle(sprintf('%s', analysis_name))
dev.off()

pca = plotPCA(vst, "Treatment", returnData=TRUE)

pdf(sprintf('%s_PCAlabeled.pdf', analysis_name))
p1 = pca %>%
  ggplot(aes(x=PC1, y=PC2, color=group)) +
  geom_point() +
  geom_text_repel(aes(label=gsub("_1.*","",rownames(pca)))) +
  ggtitle(sprintf('%s', analysis_name))
dev.off()

#We will chop off the version number of the gene IDs, so that we can better look up their annotation information later.
#However, we have a few genes which would have duplicated gene IDs after chopping off the version number, so in order to proceed we have to also use make.unique to indicate that some genes are duplicated. (It might be worth looking into why we have multiple versions of genes with the same base ID coming from our annotation.)
head(dds)
table(duplicated(substr(rownames(dds),1,18)))
rownames(dds) <- make.unique(substr(rownames(dds),1,18))
head(dds)

#Now we can run our differential expression pipeline. First, it is sometimes convenient to remove genes where all the samples have very small counts. It's less of an issue for the statistical methods, and mostly just wasted computation, as it is not possible for these genes to exhibit statistical significance for differential expression. Here we count how many genes (out of those with at least a single count) have 3 samples with a count of 10 or more:
dds <- dds[rowSums(counts(dds)) > 0,]
keep_dds <- rowSums(counts(dds) >= 1) >= 3
table(keep_dds)
dds_over1 <- dds[keep_dds,] #filter them out

dds_over1 <- DESeq(dds_over1)
resultsNames(dds_over1)
ResName <- resultsNames(dds_over1)
ResName
ResName_input <- ResName[2]
ResName_input
res_dds_over1 <- results(dds_over1, name = ResName_input)
head(res_dds_over1)
res_dds_over1.sort <- res_dds_over1[order(res_dds_over1$pvalue),]

#plotMA(res_dds_over1, ylim=c(-5,5))

summary(res_dds_over1)

#add gene names (symbols) to deseq results file
#modify your deseq results (res) table to take off numbers after decimal point to allow for matching to this database, needed to install org.Mm.eg.db previously for this to work
geneIDs <- substr(rownames(res_dds_over1), 1, 18)
# running mapIDs: collect gene symbols for the ensembl names in your geneID list
gene_symbols <- mapIds(org.Mm.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
#add gene symbols as a new column to your res file
res_dds_over1$GeneSymbol <- gene_symbols

#make volcano plot
pdf(sprintf('%s_Volcano.pdf', analysis_name))
with(res_dds_over1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))
# Add colored points: red if padj<0.05. (Other options are for orange of log2FC>1, green if both)
with(subset(res_dds_over1, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
#with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
res_dds_over1 <- as.data.frame(res_dds_over1)
res_dds_over1 <- res_dds_over1[ ,c(7, 1:6)]

write.csv((res_dds_over1),
          file=sprintf("%s_DESeqRes.csv", analysis_name))

# Generate DEG lists
signif = res_dds_over1 %>%
  filter(padj < 0.05)

down = signif %>%
  filter(log2FoldChange < 0, 
         GeneSymbol != 'NA')

up = signif %>%
  filter(log2FoldChange > 0,
         GeneSymbol != 'NA')

write_tsv(as_tibble(rownames(down)), file = sprintf("%s_down_ensembl.txt", analysis_name),
          col_names = F)
write_tsv(as_tibble(rownames(up)), file = sprintf("%s_up_ensembl.txt", analysis_name),
          col_names = F)
write_tsv(as_tibble(down$GeneSymbol), file = sprintf("%s_down_geneSymbol.txt", analysis_name),
          col_names = F)
write_tsv(as_tibble(up$GeneSymbol), file = sprintf("%s_up_geneSymbol.txt", analysis_name),
          col_names = F)

res_dds_over1 %>%
  ggplot(aes(x=baseMean)) + 
  geom_histogram(binwidth = 1) +
  xlim(0,100)

# Generate background list
background = res_dds_over1 %>%
  filter(baseMean > 5,
         GeneSymbol != 'NA') 

write_tsv(as_tibble(rownames(background)), file = sprintf("%s_background_ensembl.txt", analysis_name),
          col_names = F)
write_tsv(as_tibble(background$GeneSymbol), file = sprintf("%s_background_geneSymbol.txt", analysis_name),
          col_names = F)

