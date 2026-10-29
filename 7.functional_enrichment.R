library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
args <- commandArgs(FALSE); root <- dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)])))
out <- file.path(root,'analysis_outputs/enrichment'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
prep <- readRDS(file.path(root,'analysis_outputs/prepared/prepared.rds'))
x <- prep$expression; background <- rownames(x)[head(order(apply(x,1,var),decreasing=TRUE),5000)]
map <- function(s) unique(na.omit(mapIds(org.Hs.eg.db,keys=unique(s),keytype='SYMBOL',column='ENTREZID',multiVals='first')))
sets <- list(candidates=read.csv(file.path(root,'analysis_outputs/candidates/ranked_candidates.csv'))$gene,
 representatives=read.csv(file.path(root,'analysis_outputs/candidates/representative_subset.csv'))$gene)
de <- read.csv(file.path(root,'analysis_outputs/differential/median_group_all_genes.csv'))
sets$DE_up <- de$gene[de$significant & de$logFC>0]; sets$DE_down <- de$gene[de$significant & de$logFC<0]
for(label in names(sets)) {
 universe <- map(if(grepl('^DE_',label)) rownames(x) else background)
 query <- intersect(map(sets[[label]]),universe)
 if(length(query)<3) next
 analyses <- list(GO=enrichGO(query,OrgDb=org.Hs.eg.db,universe=universe,ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1),
 KEGG=enrichKEGG(query,organism='hsa',universe=universe,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1))
 for(kind in names(analyses)) {
  result <- analyses[[kind]]; tab <- as.data.frame(result)
  if(!nrow(tab)) next
  ratio <- function(s) vapply(strsplit(s,'/',fixed=TRUE),function(v) as.numeric(v[1])/as.numeric(v[2]),numeric(1))
  tab$FoldEnrichment <- ratio(tab$GeneRatio)/ratio(tab$BgRatio)
  write.csv(tab,file.path(out,paste0(label,'_',kind,'_all.csv')),row.names=FALSE)
  filtered <- tab[tab$p.adjust<.05 & tab$FoldEnrichment>=1.5 & tab$Count>=3,,drop=FALSE]
  write.csv(filtered,file.path(out,paste0(label,'_',kind,'_filtered.csv')),row.names=FALSE)
  if(nrow(filtered)) {
   result@result <- result@result[result@result$ID %in% filtered$ID,,drop=FALSE]
   pdf(file.path(out,paste0(label,'_',kind,'.pdf')),width=10,height=8)
   print(dotplot(result,showCategory=20)); dev.off()
  }
 }
}
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
