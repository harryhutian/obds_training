library(DESeq2)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(biomaRt)
library(EnsDb.Mmusculus.v79)

sample_table <- read_tsv("DESeq2/obds_sampletable.tsv")
counts_table <- read_tsv("DESeq2/obds_countstable.tsv.gz")

counts_table <- column_to_rownames(counts_table, "Geneid")
counts_table <- as.matrix(counts_table)

sample_table <- column_to_rownames(sample_table,"Sample_accession")

# do the colnames of counts_table  ==  to rownames of sample_table 
table(colnames(counts_table) == rownames(sample_table))

# set Egr2/3 DKO CD8 cells as the reference level
# seperate up sample_titile column

sample_table <- sample_table %>%
    separate(sample_title, c("egr_locus", "genotype", "cell_type", "replicate"), sep = "_") %>%
    unite(col = "condition", egr_locus, genotype, cell_type, sep = "_") %>%
    dplyr::select(-c(species, library_layout)) %>%
    mutate(condition = factor(condition, levels = c("Egr2/3_DKO_CD8", "Egr2/3_DKO_CD4", "Egr2_Kin_CD4", "Egr2_Kin_CD8")))

levels(sample_table$condition)
# Generate a DESeqDataSet object named dds

dds <- DESeqDataSetFromMatrix(counts_table, 
                       sample_table, 
                       ~condition
                       )

colData(dds)
rowRanges(dds) # access to the info about the row(sample info)

design(dds) # access the design formula
counts(dds) # access the counts matrix, and it is equivalent to assays(dds)$counts


# Calculate the size factors for each sample – estimateSizeFactors()

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Generate a bar plot of the size factors for each sample, coloured by condition/group


sizefactorsdf <- data.frame(sample = names(sizeFactors(dds)),
                           size_factor = sizeFactors(dds),
                           sample_group = colData(dds)$condition)

sizefactorsdf %>%
    ggplot(aes(y = size_factor, 
               x = sample, 
               fill = sample_group))+
    geom_col()+
    theme(axis.text.x=element_text(angle = 45, hjust=1))

# Obtain dispersion estimates for each gene – estimateDispersions()
dds <- estimateDispersions(dds)
dispersions(dds)

# Plot the per-gene dispersion estimates (DESeq2 has a helper function for this)

plotDispEsts(dds)

# Perform the Wald test – nbinomWaldTest()
dds <- nbinomWaldTest(dds)

# Use the DESeq() function to perform (estimatesizeFunction,estimateDispersion, Waldtest)


# Access the coefficients of the NB GLM, NAs may represent independently filtered out genes, or several other reasons

head(coef(dds),3)

# Access the results table for the comparison between CD8+ and CD4+ T cells from Egr2/3 DKO mice
res <- results(dds, contrast = c('condition', 'Egr2/3_DKO_CD4','Egr2/3_DKO_CD8'))


# Plot a histogram of the raw and BH-adjusted p-values?

res_df <- as.data.frame(res)

ggplot1 <- ggplot(res_df, aes(x= padj))+
    geom_histogram()

ggplot2 <- ggplot(res_df, aes(x = pvalue))+
    geom_histogram()

plot_grid(ggplot1, ggplot2)

# Generate an MA plot of the log2 FC values for all genes

plotMA(res) # comparing CD4 to CD8 DKO

# Shrink the log2 FC values using the normal, apeglm and ashr methods

lfcShrink(res, coef = 'condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8')

resNorm <- lfcShrink(dds, contrast = c('condition','Egr2/3_DKO_CD4','Egr2/3_DKO_CD8'), type = 'normal')

resAsh <- lfcShrink(dds, contrast = c('condition','Egr2/3_DKO_CD4','Egr2/3_DKO_CD8'), type = 'ashr')

resApeg <- lfcShrink(dds, coef = c('condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8'), type = 'apeglm')

plotMA(resAsh)
resAshplot <- recordPlot()

plotMA(resNorm)
resNormplot <- recordPlot()

plotMA(resApeg)
resApegplot <- recordPlot()

plot_grid(resAshplot, resNormplot, resApegplot)
# plotMA is base R, can store using recordplot

# Generate a results table (one shrinkage method) containing mgi symbols

resApeg_table <- as.data.frame(resApeg)

edb <- EnsDb.Mmusculus.v79
columns(edb)
gene_id <- select(edb, keys = rownames(resApeg_table), columns = c("GENEID", "SYMBOL"), keytype = 'GENEID')

# Remove all genes with a padj of NA
gene_id$GENEID[duplicated(gene_id$GENEID)]
gene_id$SYMBOL[duplicated(gene_id$SYMBOL)]

dup_genes <- dplyr::filter(gene_id, 
                    SYMBOL %in% gene_id$SYMBOL[duplicated(gene_id$SYMBOL)])


resApeg_table <- rownames_to_column(resApeg_table, "GENEID") %>%
    left_join(gene_id)

resApeg_table[is.na(resApeg_table$SYMBOL),]

sum(is.na(resApeg_table$SYMBOL))

resApeg_padj <- dplyr::filter(resApeg_table, !is.na(padj))
sum(is.na(resApeg_padj$padj))


write.csv(resApeg_padj, file = "results/resApeg_padj.csv", quote = FALSE, row.names = FALSE)


filtered_table <- dplyr::filter(resApeg_padj, padj < 0.05) %>%
                  dplyr::filter(abs(log2FoldChange) > 1)
dim(filtered_table)

write.csv(resApeg_padj, file = "results/filtered_table.csv", quote = FALSE, row.names = FALSE)

vst <- vst(dds, blind = FALSE)
rlog <- rlog(dds, blind = FALSE)

assay(vst)

vsn::meanSdPlot(assay(vst))

vsn::meanSdPlot(assay(rlog))

plotPCA(vst, ntop = nrow(vst))
plotPCA(vst, ntop = 500)

top_20 <- filtered_table[order(-abs(filtered_table$log2FoldChange)),] %>%
    slice_head(n=20) %>%
    pull(GENEID)

vst_top20 <- as.data.frame(assay(vst)) %>% 
             dplyr::filter(rownames(.) %in% top_20)

pheatmap(vst_top20, scale = "row")


