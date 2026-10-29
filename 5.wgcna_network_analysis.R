library(WGCNA)
set.seed(123)
args <- commandArgs(FALSE); root <- dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)])))
prep <- readRDS(file.path(root,'analysis_outputs/prepared/prepared.rds'))
out <- file.path(root,'analysis_outputs/wgcna'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- prep$expression
x <- x[head(order(apply(x,1,var),decreasing=TRUE),5000),,drop=FALSE]
dat <- scale(t(x))
qc <- goodSamplesGenes(dat,verbose=0)
if(!qc$allOK) stop('WGCNA input failed quality checks')
soft <- pickSoftThreshold(dat,powerVector=c(1:10,seq(12,20,2)),networkType='signed',verbose=0)
write.csv(soft$fitIndices,file.path(out,'soft_threshold_diagnostics.csv'),row.names=FALSE)
power <- soft$powerEstimate
if(!is.finite(power)) stop('No supported soft power; inspect diagnostics')
net <- blockwiseModules(dat,power=power,networkType='signed',TOMType='signed',minModuleSize=50,mergeCutHeight=.25,numericLabels=FALSE,maxBlockSize=6000,verbose=0)
me <- net$MEs; sf2 <- prep$response$SF2
r <- cor(me,sf2); p <- corPvalueStudent(r,nrow(dat))
module <- data.frame(module=colnames(me),r=as.vector(r),p=as.vector(p),FDR=p.adjust(as.vector(p),'BH'))
write.csv(module,file.path(out,'module_trait_associations.csv'),row.names=FALSE)
mm <- signedKME(dat,me); gs <- as.vector(cor(dat,sf2))
genes <- data.frame(gene=colnames(dat),module=net$colors,GS=gs,kME=NA_real_)
for(i in seq_len(nrow(genes))) genes$kME[i] <- mm[i,match(paste0('kME',genes$module[i]),colnames(mm))]
genes$screened <- genes$kME>.85 & abs(genes$GS)>.3
write.csv(genes,file.path(out,'module_membership.csv'),row.names=FALSE)
candidate <- read.csv(file.path(root,'analysis_outputs/candidates/ranked_candidates.csv'))$gene
write.csv(genes[genes$gene %in% candidate,,drop=FALSE],file.path(out,'candidate_module_membership.csv'),row.names=FALSE)
for(m in sub('^ME','',module$module[module$FDR<.05])) {
 if(m=='grey') next
 block <- dat[,net$colors==m,drop=FALSE]
 tom <- TOMsimilarityFromExpr(block,power=power,networkType='signed',TOMType='signed',verbose=0)
 exportNetworkToCytoscape(tom,edgeFile=file.path(out,paste0(m,'.edge')),nodeFile=file.path(out,paste0(m,'.node')),weighted=TRUE,threshold=.1,nodeNames=colnames(block))
}
pdf(file.path(out,'network_diagnostics.pdf'),width=10,height=7)
plot(soft$fitIndices[,1],-sign(soft$fitIndices[,3])*soft$fitIndices[,2],xlab='Soft power',ylab='Signed scale-free fit')
plotDendroAndColors(net$dendrograms[[1]],net$colors[net$blockGenes[[1]]],'Module',dendroLabels=FALSE)
barplot(module$r,names.arg=module$module,las=2,ylab='Correlation with SF2')
dev.off()
saveRDS(net,file.path(out,'network.rds')); writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
