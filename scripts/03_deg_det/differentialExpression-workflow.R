# Performing differential expression analysis with DESeq2 workflow (v1.26.0)
# Plots made using scripts from: https://gist.github.com/stephenturner/f60c1934405c127f09a6
# Ran separately on 5Y and sequin counts

library(DESeq2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(pcaExplorer)
library(ggplot2)
library(dplyr)
setwd("../output/")

countdata=read.csv('./input_data/OUT.transcript_model_grouped_counts.tsv',header = T,sep = '\t',row.names = 1)
# this file contain the normal transcripts and the novel discovered isoforms
head(countdata)
countdata <-countdata[,c(1,2,3,5,6,7)]
countdata <- as.matrix(countdata)

# Assign conditions
(condition <- factor(c(rep("Aposymbiotic", 3), rep("Symbiotic", 3))))

# Make DESeq dataset
(coldata <- data.frame(row.names=colnames(countdata), condition))
countdata<-round(countdata)
coldata$sample=rownames(coldata)
coldata$batch=c(1,2,1,2,2,2)#c(1,2,1,3,2,2,2,3)
coldata$batch=as.factor(coldata$batch)
coldata$date=as.factor(c(1,2,1,2,4,2)) #c(1,2,1,3,2,4,2,3)
coldata$ont=as.factor(c(1,1,1,1,1,1)) #c(1,1,1,2,1,1,1,2)
coldata
coldata$batch_combind=as.factor(c(1,2,1,2,3,2))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~batch+condition)
dds

# Optional filtering step to remove very low counts (chosen minimum of 5)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Run DESeq2 pipeline
dds <- DESeq(dds)
res <- DESeq2::results(dds)
res

#DESeq2 results
table(res$padj<0.05)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

resdata=resdata[resdata$padj<0.05,]
resdata=resdata[abs(resdata$log2FoldChange)>1,]
resdata=resdata[!is.na(resdata$log2FoldChange),]
resdata$log2FoldChange
table(resdata$log2FoldChange< (-1))
table(resdata$log2FoldChange> (1))

# load the transcript to gene mapping file
mapping=read.csv('../DEG_analysis/tx2gene.tsv',sep = '\t',header = F)
colnames(mapping)=c('transcript_id','Gene')
colnames(resdata)[1]='transcript_id'
resdata=merge(resdata,mapping,by.x='transcript_id',by.y='transcript_id',all.x=TRUE)
head(resdata)

mapping=read.csv('../GO_analysis/isoquant_anno/transcript_level_GO.tsv',sep = '\t')
head(mapping)
head(resdata)
mapping=mapping[,c(1,2,9,10)]
resdata=merge(resdata,mapping,by.x='transcript_id',by.y='transcript_id',all.x=TRUE)
head(resdata)

resdata=resdata[,-c(15)]
resdata$deg_result <- ifelse(resdata$log2FoldChange > 1, "up", "down")
test=resdata[duplicated(resdata$Gene),]
test=test[order(test$Gene),]
test=test[,c(1,14,15,17)]
test
write.csv(resdata, file="../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_transcript_level_raw.csv")


library(dplyr)

conflict_rows <- resdata %>%
  group_by(Gene) %>%
  filter(n() > 1, n_distinct(deg_result) > 1) %>%
  arrange(Gene, transcript_id) %>%
  ungroup()

## summary
## totally 501 transcripts, 453 genes 
## totally 10 genes shows conflict differential expression result at the transcript level
# Write DE results 
write.csv(resdata, file="../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model.csv")

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot", genecol = "black", fitcol = "cyan", finalcol = "blue", legend = TRUE)

# RLD for viewing
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Plot residual p-values
hist(res$pvalue, breaks=50, col="grey")

#Set colours for plotting
mycols <- brewer.pal(8, "Accent")[1:length(unique(condition))]

# Heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png("heatmap-samples-6samples_update.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# PCA
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomright", main="Principal Component Analysis", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1: ",as.character(pc1var),"% variance")
  pc2lab <- paste0("PC2: ",as.character(pc2var),"% variance")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
png("pca-1-2-6samples-update.png", 1000, 1000, pointsize=30)
rld_pca(rld, colors=mycols, intgroup="condition")
dev.off()

# MA Plot
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="blue", pch=20, cex=1))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot-0.05-6samples-update.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

# Volcano Plot
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, xlab="log2(Fold Change)", legendpos="topright", labelsig=FALSE, textcx=1.5, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlab=xlab, cex.axis=1.8, cex.lab=1.5, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p-adj<",sigthresh,sep=""), paste("|log2(FC)|>",lfcthresh,sep=""), "both"), cex=1.5, pch=20, col=c("blue","orange","green"))
}
pdf("diffexpr-volcanoplot-hi-res-6samples-update.pdf", 18, 15, pointsize=20)
volcanoplot(resdata, lfcthresh=2, sigthresh=0.05, xlim=c(-6, 6), ylim=c(0,33), legendpos="topright")
dev.off()

#### below is to do the gene level deg analysis #####
### map the transcript id back to the gene id level first
countdata=read.csv('./input_data/OUT.transcript_model_grouped_counts.tsv',header = T,sep = '\t',row.names = 1)
# this file contain the normal transcripts and the novel discovered isoforms
dim(countdata)
head(countdata)
countdata$transcript_id=rownames(countdata)
mapping=read.csv('../DEG_analysis/tx2gene.tsv',sep = '\t',header = F)
colnames(mapping)=c('transcript_id','Gene')
head(mapping)
countdata=merge(countdata,mapping,by='transcript_id',all.x=TRUE)
dim(countdata)
library(dplyr)

gene_countdata <- countdata %>%
  #mutate(Gene = sub("_gene$", "", Gene)) %>%   # remove only trailing "_gene"
  group_by(Gene) %>%
  summarise(across(A1:H4, sum, na.rm = TRUE), .groups = "drop")

gene_countdata=as.data.frame(gene_countdata)
rownames(gene_countdata)=gene_countdata$Gene
head(gene_countdata)
countdata=gene_countdata[,c(2,3,4,5,6,7,8)]
#remove extra columns of data from featureCounts if needed
countdata <-countdata[,c(1,2,3,5,6,7)]
countdata <- as.matrix(countdata)

# Assign conditions
(condition <- factor(c(rep("Aposymbiotic", 3), rep("Symbiotic", 3))))

# Make DESeq dataset
(coldata <- data.frame(row.names=colnames(countdata), condition))
countdata<-round(countdata)
coldata$sample=rownames(coldata)
coldata$batch=c(1,2,1,2,2,2)#c(1,2,1,3,2,2,2,3)
coldata$batch=as.factor(coldata$batch)
coldata$date=as.factor(c(1,2,1,2,4,2)) #c(1,2,1,3,2,4,2,3)
coldata$ont=as.factor(c(1,1,1,1,1,1)) #c(1,1,1,2,1,1,1,2)
coldata
coldata$batch_combind=as.factor(c(1,2,1,2,3,2))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~batch+condition)
dds

# Optional filtering step to remove very low counts (chosen minimum of 5)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Run DESeq2 pipeline
dds <- DESeq(dds)
res <- DESeq2::results(dds)
res


#DESeq2 results
table(res$padj<0.05)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_gene_level_raw.csv")

resdata=resdata[resdata$padj<0.05,]
dim(resdata)
resdata=resdata[abs(resdata$log2FoldChange)>1,]
resdata=resdata[!is.na(resdata$log2FoldChange),]
resdata$log2FoldChange
table(resdata$log2FoldChange< (-1))
table(resdata$log2FoldChange> (1))

## summary
## totally 456 genes, 245 up regulated, 211 down regulated 
## 39 are novol genes, 8.6%

# Write DE results 
write.csv(resdata, file="../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_gene_level.csv")

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot", genecol = "black", fitcol = "cyan", finalcol = "blue", legend = TRUE)

# RLD for viewing
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Plot residual p-values
hist(res$pvalue, breaks=50, col="grey")

#Set colours for plotting
mycols <- brewer.pal(8, "Accent")[1:length(unique(condition))]

# Heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png("heatmap-samples-6samples_update_gene_level.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# PCA
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomright", main="Principal Component Analysis", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1: ",as.character(pc1var),"% variance")
  pc2lab <- paste0("PC2: ",as.character(pc2var),"% variance")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
png("pca-1-2-6samples-update-gene_level.png", 1000, 1000, pointsize=30)
rld_pca(rld, colors=mycols, intgroup="condition")
dev.off()

# MA Plot
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="blue", pch=20, cex=1))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot-0.05-6samples-update-gene_level.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

# Volcano Plot
volcanoplot_pub <- function(
    res,
    lfcthresh = 1,
    sigthresh = 0.05,
    xlab = expression(log[2]("Fold Change")),
    ylab = expression(-log[10]("adjusted p-value")),
    main = NULL,
    legendpos = "topright",
    labelsig = FALSE,
    label_col = "Gene",
    max_labels = 15,
    point_cex = 1.0,
    axis_cex = 1.4,
    lab_cex = 1.6,
    main_cex = 1.6
) {
  # keep only valid rows
  df <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  # avoid Inf when padj = 0
  df$padj_plot <- pmax(df$padj, .Machine$double.xmin)
  df$neglog10_padj <- -log10(df$padj_plot)
  
  # classify points
  df$group <- "NS"
  df$group[df$padj < sigthresh] <- "padj"
  df$group[abs(df$log2FoldChange) > lfcthresh] <- "LFC"
  df$group[df$padj < sigthresh & abs(df$log2FoldChange) > lfcthresh] <- "Both"
  
  # colors
  cols <- c(
    "NS"   = "grey75",
    "padj" = "#3B82F6",   # blue
    "LFC"  = "#F59E0B",   # orange
    "Both" = "#DC2626"    # red
  )
  
  # plotting limits
  xlim <- range(df$log2FoldChange, na.rm = TRUE)
  ylim <- c(0, max(df$neglog10_padj, na.rm = TRUE) * 1.05)
  
  # draw plot
  plot(
    df$log2FoldChange, df$neglog10_padj,
    pch = 16,
    col = cols[df$group],
    cex = point_cex,
    xlab = xlab,
    ylab = ylab,
    main = main,
    xlim = xlim,
    ylim = ylim,
    cex.axis = axis_cex,
    cex.lab = lab_cex,
    cex.main = main_cex,
    las = 1
  )
  
  # threshold lines
  abline(v = c(-lfcthresh, lfcthresh), lty = 2, lwd = 1.5, col = "grey40")
  abline(h = -log10(sigthresh), lty = 2, lwd = 1.5, col = "grey40")
  
  # legend
  legend(
    legendpos,
    legend = c(
      "Not significant",
      paste0("padj < ", sigthresh),
      paste0("|log2FC| > ", lfcthresh),
      "Both"
    ),
    col = cols[c("NS", "padj", "LFC", "Both")],
    pch = 16,
    pt.cex = 1.2,
    bty = "n",
    cex = 1.2
  )
  
  # optional labels for top significant genes
  if (labelsig && label_col %in% colnames(df)) {
    lab_df <- df[df$group == "Both", ]
    lab_df <- lab_df[order(lab_df$padj, -abs(lab_df$log2FoldChange)), ]
    lab_df <- head(lab_df, max_labels)
    
    with(
      lab_df,
      text(
        x = log2FoldChange,
        y = neglog10_padj,
        labels = lab_df[[label_col]],
        pos = 3,
        cex = 0.9,
        xpd = NA
      )
    )
  }
}

png(
  filename = "diffexpr_volcanoplot_gene_level.png",
  width = 7,
  height = 6,
  units = "in",
  res = 600
)

par(
  mar = c(5.5, 5.5, 2, 1),
  mgp = c(3, 1, 0),
  family = "sans"
)

volcanoplot_pub(
  resdata,
  lfcthresh = 1,
  sigthresh = 0.05,
  main = "Differential expression",
  labelsig = FALSE,
  label_col = "Gene"
)

dev.off()

#########################################################################
##### now comparing the deg results and the de transcript results #######
#########################################################################
deg=read.csv('../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_gene_level_raw.csv',header = T,row.names = 1)
det=read.csv('../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_transcript_level_raw.csv',header = T,row.names = 1)

deg=deg[deg$padj<0.05,]
det=det[det$padj<0.05,]

deg=deg[abs(deg$log2FoldChange)>1,]
det=det[abs(det$log2FoldChange)>1,]
deg=deg[!is.na(deg$log2FoldChange),]
dim(deg)

det=det[!is.na(det$log2FoldChange),]
dim(det)

head(deg)
head(det)
deg$deg_result <- ifelse(deg$log2FoldChange > 1, "up", "down")
library(VennDiagram)
library(grid)

genes_deg <- unique(na.omit(deg$Gene))
genes_det <- unique(na.omit(det$Gene))

venn.plot <- venn.diagram(
  x = list(
    DEG = genes_deg,
    DET = genes_det
  ),
  filename = NULL,
  fill = c("#4C78A8", "#F58518"),
  alpha = 0.5,
  cex = 1.8,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

png("venn_deg_det_genes.png", width = 6, height = 6, units = "in", res = 600)
grid.newpage()
grid.draw(venn.plot)
dev.off()


library(UpSetR)
gene_list <- list(
  DEG = genes_deg,
  DET = genes_det
)

png("upset_deg_det_genes.png", width = 7, height = 5, units = "in", res = 600)
upset(fromList(gene_list), order.by = "freq")
dev.off()


library(dplyr)

tx_gene_summary <- det %>%
  filter(!is.na(Gene), !is.na(deg_result)) %>%
  group_by(Gene) %>%
  dplyr::summarise(
    tx_direction = case_when(
      all(deg_result == "up") ~ "up",
      all(deg_result == "down") ~ "down",
      TRUE ~ "mixed"
    ),
    n_transcripts = n(),
    .groups = "drop"
  )

tx_per_gene <- det %>%
  count(Gene, name = "n_transcripts")

dist_df <- tx_per_gene %>%
  count(n_transcripts, name = "gene_count") %>%
  mutate(
    percent = gene_count / sum(gene_count) * 100,
    label = paste0(gene_count, " (", round(percent, 1), "%)")
  )

dist_df

p <- ggplot(dist_df, aes(x = factor(n_transcripts), y = gene_count)) +
  geom_col(width = 0.7, fill = "steelblue") +
  geom_text(aes(label = label), vjust = -0.3, size = 7) +
  labs(
    x = "Number of DET transcripts per gene",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 24) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  expand_limits(y = max(dist_df$gene_count) * 1.12)
p

ggsave(
  "det_transcripts_per_gene_barplot.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 600
)

tx_per_gene_dir <- det %>%
  group_by(Gene) %>%
  summarise(
    n_transcripts = n(),
    direction = case_when(
      all(deg_result == "up") ~ "up",
      all(deg_result == "down") ~ "down",
      TRUE ~ "mixed"
    ),
    .groups = "drop"
  )

dist_dir <- tx_per_gene_dir %>%
  count(n_transcripts, direction)

p3 <- ggplot(dist_dir, aes(x = factor(n_transcripts), y = n, fill = direction)) +
  geom_col(position = "stack", width = 0.7) +
  labs(
    x = "Number of DET transcripts per gene",
    y = "Number of genes",
    fill = "Direction"
  ) +
  theme_classic(base_size = 24) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

p3
ggsave(
  "det_transcripts_per_gene_with_direction_barplot.png",
  plot = p3,
  width = 7,
  height = 5,
  dpi = 600
)

gene_summary <- deg %>%
  filter(!is.na(Gene), !is.na(deg_result)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  dplyr::select(Gene, gene_direction = deg_result)

compare_df <- full_join(gene_summary, tx_gene_summary, by = "Gene")
head(compare_df)
table(compare_df$gene_direction, compare_df$tx_direction, useNA = "ifany")

library(ggplot2)

plot_df <- compare_df %>%
  filter(!is.na(gene_direction), !is.na(tx_direction)) %>%
  count(gene_direction, tx_direction)

png("deg_det_direction_heatmap.png", width = 6, height = 5, units = "in", res = 600)
ggplot(plot_df, aes(x = tx_direction, y = gene_direction, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 5) +
  scale_fill_gradient(low = "grey90", high = "steelblue") +
  labs(
    x = "Transcript-level gene summary",
    y = "Gene-level DEG direction",
    title = "Agreement between DEG and DET directions"
  ) +
  theme_classic(base_size = 14)
dev.off()


plot_df2 <- compare_df %>%
  filter(!is.na(gene_direction), !is.na(tx_direction)) %>%
  count(gene_direction, tx_direction) %>%
  group_by(gene_direction) %>%
  mutate(prop = n / sum(n))

png("deg_det_direction_barplot.png", width = 6, height = 5, units = "in", res = 600)
ggplot(plot_df2, aes(x = gene_direction, y = prop, fill = tx_direction)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Gene-level DEG direction",
    y = "Proportion",
    fill = "Transcript summary",
    title = "Transcript-level direction within DEG genes"
  ) +
  theme_classic(base_size = 14)
dev.off()

summary_table <- compare_df %>%
  mutate(
    comparison = case_when(
      is.na(gene_direction) & !is.na(tx_direction) ~ "DET only",
      !is.na(gene_direction) & is.na(tx_direction) ~ "DEG only",
      gene_direction == tx_direction ~ "Same direction",
      tx_direction == "mixed" ~ "Mixed transcripts",
      TRUE ~ "Opposite direction"
    )
  )

table(summary_table$comparison)

### now get the unique genes from both file
library(dplyr)
deg_only <- deg %>%
  anti_join(det %>% distinct(Gene), by = "Gene")

det_only <- det %>%
  anti_join(deg %>% distinct(Gene), by = "Gene")

det_only_summary=det_only %>%
  filter(!is.na(Gene), !is.na(deg_result)) %>%
  group_by(Gene) %>%
  dplyr::summarise(
    tx_direction = case_when(
      all(deg_result == "up") ~ "up",
      all(deg_result == "down") ~ "down",
      TRUE ~ "mixed"
    ),
    n_transcripts = n(),
    .groups = "drop"
  )
table(det_only_summary$tx_direction)

write.csv(deg_only$Gene,'unique_genes_detected_by_deg_analysis.txt',row.names = F,col.names =NA,quote = F)
write.csv(det_only$Gene,'unique_genes_detected_by_det_analysis.txt',row.names = F,col.names =NA,quote = F)

det$transcript_id
write.csv(det[det$deg_result=='up',]$transcript_id,'det_up_regulated_sym_analysis.txt',row.names = F,col.names =NA,quote = F)
write.csv(det[det$deg_result=='down',]$transcript_id,'det_down_regulated_sym_analysis.txt',row.names = F,quote = F)
write.csv(det_only_summary[det_only_summary$tx_direction=='mixed',]$Gene, 'unique_genes_from_det_with_mixed_directions_transcripts.txt',row.names = F,quote = F)
