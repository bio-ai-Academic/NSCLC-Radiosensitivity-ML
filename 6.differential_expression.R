library(limma)
set.seed(123)
args <- commandArgs(FALSE); root <- dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)])))
prep <- readRDS(file.path(root,'analysis_outputs/prepared/prepared.rds'))
out <- file.path(root,'analysis_outputs/differential'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
candidate <- read.csv(file.path(root,'analysis_outputs/candidates/ranked_candidates.csv'))$gene
results <- list()
for(scheme in c('median_group','extreme_group','fixed_group')) {
 labels <- prep$groups[[scheme]]; keep <- !is.na(labels)
 labels <- sub('Highly ','',labels[keep]); group <- factor(labels,levels=c('Sensitive','Resistant'))
 if(any(table(group)<3)) stop('Insufficient replicates in a group')
 meta <- prep$metadata[keep,,drop=FALSE]; covariates <- data.frame(group=group)
 for(v in intersect(c('histology','batch'),names(meta))) {
  if(anyNA(meta[[v]])) stop(paste('Missing covariate',v))
  if(length(unique(meta[[v]]))>1) covariates[[v]] <- factor(meta[[v]])
 }
 design <- model.matrix(reformulate(names(covariates)),covariates)
 if(qr(design)$rank<ncol(design)) stop('Group and covariates confounded; inspect design')
 fit <- eBayes(lmFit(prep$expression[,keep,drop=FALSE],design),trend=TRUE)
 tab <- topTable(fit,coef='groupResistant',number=Inf,sort.by='none')
 tab$gene <- rownames(tab); tab$significant <- abs(tab$logFC)>=1 & tab$adj.P.Val<.05
 results[[scheme]] <- tab
 write.csv(tab,file.path(out,paste0(scheme,'_all_genes.csv')),row.names=FALSE)
 write.csv(tab[tab$significant,,drop=FALSE],file.path(out,paste0(scheme,'_significant.csv')),row.names=FALSE)
 write.csv(as.data.frame(table(group)),file.path(out,paste0(scheme,'_counts.csv')),row.names=FALSE)
 pdf(file.path(out,paste0(scheme,'_plots.pdf')),width=9,height=6)
 plot(tab$logFC,-log10(pmax(tab$adj.P.Val,.Machine$double.xmin)),col=ifelse(tab$significant,'red','grey'),pch=16,cex=.5,xlab='Log2 fold change: Resistant / Sensitive',ylab='-log10 FDR')
 plotMA(fit,coef='groupResistant')
 top <- head(tab$gene[order(tab$adj.P.Val)],50)
 heatmap(prep$expression[top,keep,drop=FALSE],scale='row',labCol=labels,labRow=NA)
 dev.off()
}
a <- results$median_group; b <- results$extreme_group
sa <- a$gene[a$significant]; sb <- b$gene[b$significant]
write.csv(data.frame(jaccard=if(length(union(sa,sb))) length(intersect(sa,sb))/length(union(sa,sb)) else NA_real_,logFC_spearman=cor(a$logFC,b$logFC,method='spearman')),file.path(out,'threshold_consistency.csv'),row.names=FALSE)
write.csv(data.frame(candidate_only=length(setdiff(candidate,sa)),shared=length(intersect(candidate,sa)),DE_only=length(setdiff(sa,candidate))),file.path(out,'candidate_overlap_counts.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
