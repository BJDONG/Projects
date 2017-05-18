# MyScRNA
library(vioplot)
library(scater)
library(scran)
library(ggplot2)
library(ComplexHeatmap)
library(amap)

initialization <- function(pheno){
  color <- c()
  for (i in 1:length(pheno$Identity)){
    if (pheno$Pro[i] == 'A') color[i] = 'red'
    if (pheno$Pro[i] == 'B') color[i] = 'deepskyblue'
    if (pheno$Pro[i] == 'C') color[i] = 'darkmagenta'
    if (pheno$Pro[i] == 'D') color[i] = 'orange'
    if (pheno$Pro[i] == 'E') color[i] = 'seagreen'
    if (pheno$Pro[i] == 'F') color[i] = 'grey'
  }
}

findHighQualitySamples <- function(qc_data, lbu = 300000, lbl = 2000000, 
                                   mr = 30, ngu = 3000, ngl = 11000 ){
  clean_qc <- qc_data[which(qc_data$effective_reads > lbu & qc_data$effective_reads < lbl
                            & qc_data$mapping_rate > mr
                            & qc_data$genes > ngu & qc_data$genes < ngl),]
  par(mfrow=c(2,2), cex.lab = 1.5, pin = c(4, 2))
  hist(qc$total_reads/1e6, xlab = 'Clean reads (millions)', main = '', 
       breaks = 20, col = 'grey80', ylab = 'Number of cells')
  hist(qc$effective_reads/1e6, xlab = 'Library size (millions)', main = '', 
       breaks = 20, col = 'grey80', ylab = 'Number of cells')
  abline(v = c(lbu/1e6, lbl/1e6), lty = 2, col = 'red')
  hist(qc$mapping_rate, xlab = 'Mapped rates', main = '', 
       breaks = 20, col = 'grey80', ylab = 'Number of cells')
  abline(v = mr, lty = 2, col = 'red')
  hist(qc$genes, xlab = 'Number of expressed genes', main = '', 
       breaks = 20, col = 'grey80', ylab = 'Number of cells')
  abline(v = c(ngu, ngl), lty = 2, col = 'red')
  samples <- rownames(clean_qc)
  return(samples)
}

findHighQualityGenes <- function(umi){
  sce <- newSCESet(countData=umi)
  numcells <- nexprs(sce, byrow=TRUE)
  # genes expressed in more than 10% of sequenced samples
  n <- as.integer(length(umi)/10)
  alt.keep <- numcells >= n
  genes <- rownames(umi)[alt.keep]
  return(genes)
}

findHVG <- function(umi){
  sce <- newSCESet(countData=all.counts)
  sce <- computeSumFactors(sce, sizes=c(20, 40, 60, 80))
  summary(sizeFactors(sce))
  sce <- normalize(sce)
  var.fit <- trendVar(sce, trend="loess", use.spikes=FALSE, span=0.2)
  var.out <- decomposeVar(sce, var.fit)
  hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
  hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
}

pca <- function(tpm, pheno, color = color){
  data.pca <- prcomp(tpm)
  pc1 <- data.pca$x[ ,1]
  pc2 <- data.pca$x[ ,2]
  par(fig=c(0, 0.8, 0, 1))
  plot(pc1, pc2, col = color, pch = 19, cex.lab = 1.3,
       xlab = 'PC1:15.46% variance', ylab = 'PC2:5.18% variance')
  par(fig = c(0.7, 1, 0, 1), new=TRUE)
  legend('top', title = 'Clusters', 
         legend = c('A', 'B', 'C', 'D', 'E', 'F'), 
         pch = c(19, 19, 19, 19, 19, 19), 
         box.col = 'white', text.font = 1,
         col = c('red', 'deepskyblue', 'darkmagenta', 'orange', 'seagreen', 'grey'))
}

tsne <- function(tpm, pheno, color = color){
  par(fig=c(0, 0.8, 0, 1))
  plot(tsne1, tsne2, col = color, pch = 19, cex.lab = 1.3,
       xlab = 'Dimension1', ylab = 'Dimension2')
  par(fig = c(0.7, 1, 0, 1), new=TRUE)
  legend('top', title = 'Clusters', 
         legend = c('A', 'B', 'C', 'D', 'E', 'F'),
         pch = c(19, 19, 19, 19, 19, 19), 
         box.col = 'white', text.font = 1,
         col = c('red', 'deepskyblue', 'darkmagenta', 'orange', 'seagreen', 'grey'))
}

pseudoTime <- function(umi, pheno, gene){
  gene_short_name <- row.names(umi)
  feature <- as.data.frame(gene_short_name)
  rownames(feature) <- gene_short_name
  pd <- new("AnnotatedDataFrame", data = pheno)
  fd <- new("AnnotatedDataFrame", data = feature)
  MONO <- newCellDataSet(as.matrix(expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=negbinomial())
  MONO <- estimateSizeFactors(MONO)
  MONO <- estimateDispersions(MONO)
  disp_table <- dispersionTable(MONO)
  MONO <- setOrderingFilter(MONO, gene)
  set.seed(1234)
  MONO <- reduceDimension(MONO, max_components=2)
  MONO <- orderCells(MONO, reverse=TRUE)
  plot_cell_trajectory(MONO, color_by="State", cell_size = 2, show_branch_points = FALSE)
}

plotVio <- function(tpm, pheno, gene, outpath = getwd(), type = 'mouse'){
  temp1 <- as.data.frame(t(tpm[pheno[which(pheno$Pro == 1), 1]]))
  temp2 <- as.data.frame(t(tpm[pheno[which(pheno$Pro == 2), 1]]))
  temp3 <- as.data.frame(t(tpm[pheno[which(pheno$Pro == 3), 1]]))
  temp4 <- as.data.frame(t(tpm[pheno[which(pheno$Pro == 4), 1]]))
  temp5 <- as.data.frame(t(tpm[pheno[which(pheno$Pro == 5 | pheno$Pro == 6), 1]]))
  cluster1 <- as.matrix(subset(temp1, ,colnames(temp1) == gene))
  cluster2 <- as.matrix(subset(temp2, ,colnames(temp2) == gene))
  cluster3 <- as.matrix(subset(temp3, ,colnames(temp3) == gene))
  cluster4 <- as.matrix(subset(temp4, ,colnames(temp4) == gene))
  cluster5 <- as.matrix(subset(temp5, ,colnames(temp5) == gene))
  if (sum(cluster1) != 0 & sum(cluster2)!=0 & sum(cluster3)!=0 & sum(cluster4)!=0 & sum(cluster5)!=0){
    png(filename = paste(output, '/', gene, '.png', sep = ''), width = 600, height = 400, type = 'cairo')
    vioplot(log10(cluster1+1), log10(cluster2+1), log10(cluster3+1), log10(cluster4+1),
            log10(cluster5+1), 
            names = c('cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5'),
            col = 'grey')
    title(main = gene, ylab = 'Expression level (log10(TPM))', cex.main = 1.5, cex.lab = 1.2)
    dev.off()
  }
}
