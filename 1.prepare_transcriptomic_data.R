set.seed(123)
args <- commandArgs(FALSE)
root <- dirname(normalizePath(sub('^--file=', '', args[grepl('^--file=', args)])))
input <- file.path(root, 'analysis_inputs')
out <- file.path(root, 'analysis_outputs', 'prepared')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
expr <- read.csv(file.path(input,'expression_log2.csv'),row.names=1,check.names=FALSE)
response <- read.csv(file.path(input,'radiation_response.csv'),check.names=FALSE)
meta <- read.csv(file.path(input,'sample_metadata.csv'),check.names=FALSE)
stopifnot(all(c('sample_id','SF2','AUC_recomputed') %in% names(response)))
stopifnot(all(c('sample_id','histology') %in% names(meta)))
if(anyDuplicated(meta$sample_id)) stop('Duplicate metadata sample IDs')
response$sample_id <- sub('_radiation_[0-9]+$', '', response$sample_id)
metrics <- intersect(c('SF2','AUC_published','AUC_recomputed','alpha','beta'),names(response))
response[metrics] <- lapply(response[metrics],function(v) suppressWarnings(as.numeric(v)))
response$SF2[response$SF2<0 | response$SF2>1] <- NA_real_
for(v in intersect(c('AUC_published','AUC_recomputed'),metrics)) response[[v]][response[[v]]<=0] <- NA_real_
response <- aggregate(response[metrics],list(sample_id=response$sample_id),function(v) if(all(is.na(v))) NA_real_ else median(v,na.rm=TRUE))
ids <- Reduce(intersect,list(colnames(expr),response$sample_id,meta$sample_id))
if(length(ids)<15) stop('Fewer than 15 matched samples; check input IDs')
allids <- Reduce(union,list(colnames(expr),response$sample_id,meta$sample_id))
write.csv(data.frame(sample_id=allids,expression=allids %in% colnames(expr),response=allids %in% response$sample_id,metadata=allids %in% meta$sample_id,matched=allids %in% ids),file.path(out,'matching_log.csv'),row.names=FALSE)
expr <- as.matrix(expr[,ids,drop=FALSE]); storage.mode(expr) <- 'numeric'
if(any(!is.finite(expr))) stop('Expression must contain finite normalized log2 values')
response <- response[match(ids,response$sample_id),,drop=FALSE]
meta <- meta[match(ids,meta$sample_id),,drop=FALSE]
if(any(!is.finite(response$SF2))) stop('Observed SF2 outcomes required; targets are not imputed')
keep <- rowMeans(expr>0)>=.1 & apply(expr,1,var)>0
expr <- expr[keep,,drop=FALSE]
thresholds <- quantile(response$SF2,c(.25,.5,.75),names=FALSE)
groups <- data.frame(sample_id=ids,SF2=response$SF2,
 median_group=ifelse(response$SF2<=thresholds[2],'Sensitive','Resistant'),
 extreme_group=ifelse(response$SF2<=thresholds[1],'Highly Sensitive',ifelse(response$SF2>=thresholds[3],'Highly Resistant',NA)),
 fixed_group=ifelse(response$SF2<=.5,'Sensitive','Resistant'))
write.csv(expr,file.path(out,'expression_log2.csv'))
write.csv(response,file.path(out,'response.csv'),row.names=FALSE)
write.csv(meta,file.path(out,'metadata.csv'),row.names=FALSE)
write.csv(groups,file.path(out,'groups.csv'),row.names=FALSE)
write.csv(data.frame(threshold=c('Q1','Median','Q3'),value=thresholds),file.path(out,'thresholds.csv'),row.names=FALSE)
write.csv(as.data.frame(table(groups$median_group)),file.path(out,'group_counts.csv'),row.names=FALSE)
pdf(file.path(out,'quality_control.pdf'),width=9,height=7)
boxplot(expr,outline=FALSE,las=2,cex.axis=.5,ylab='Log2 expression')
pca <- prcomp(t(expr),scale.=FALSE)
plot(pca$x[,1:2],col=as.integer(factor(groups$median_group)),pch=19,xlab='PC1',ylab='PC2')
hist(response$SF2,main='Observed SF2',xlab='SF2')
dev.off()
saveRDS(list(expression=expr,response=response,metadata=meta,groups=groups),file.path(out,'prepared.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
