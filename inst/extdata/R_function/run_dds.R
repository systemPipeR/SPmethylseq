####################################################################
## Run DSS with cytosine methylation report from Bismark ##
####################################################################
run_dss <- function(targets, args, cmp=cmp[[1]], smoothing=TRUE) {
	##Settings
  if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) 
  samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
  targets_list <- as.list(targetsout(args)$FileName)
  names(targets_list) <- names(samples)	
  dss_input <- list() 
  
  for (i in seq(along=targets_list)){
    file.list <- as.character(list.files(paste0(targets_list[i]), "*.cov", full.names=TRUE))
    CpG_report <- read.table(file.list, col.names=c("chr","start","end","M_percentage", "countM","countUnM"))
    CpG_report <- data.frame(chr=CpG_report$chr, pos=CpG_report$start, N=CpG_report$countM+CpG_report$countUnM, X=CpG_report$countM) 
    dss_input <- c(dss_input, list(CpG_report))
  }
  names(dss_input) <- names(samples) 
  
  dmlDF <- list()
  for (i in seq(along=cmp[,1])) {
    group1 <- samples[samples %in% cmp[i,][1]]
    group2 <- samples[samples %in% cmp[i,][2]]
    sub <- c(dss_input[names(dss_input) %in% names(samples[samples %in% cmp[i,][1]])], dss_input[names(dss_input) %in% names(samples[samples %in% cmp[i,][2]])])   
    ##Create BSseq class
    dss.sub.samples <- c(names(group1), names(group2))
    dss.obj <- DSS::makeBSseqData(sub, dss.sub.samples)
    saveRDS(dss.obj, file=paste0("./results/dss.input_", paste0(cmp[i,][1], sep="-", cmp[i,][2]), ".rds"))
    ##Differntially methylated loci (DML) for two group comparisons
    if(smoothing==TRUE){
      dmlTest <- DSS::DMLtest(dss_input, group1=names(group1), group2=names(group2), smoothing=TRUE)}
    if(smoothing==FALSE){
      dmlTest <- DSS::DMLtest(dss_input, group1=names(group1), group2=names(group2))}
    dmlDF <- c(dmlDF, list(dmlTest))
  }
  names <- apply(cmp, 1, paste, collapse="-")
  names(dmlDF) <- names
  return(dmlDF)
} 




## Usage:
# targetspath <- "targets.txt"
# targets <- read.delim(targetspath, comment="#")
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# args <- systemArgs(sysma="param/bismark_methyl_extractor.param",
#                    mytargets="targets_bismark.txt")
# dss.DF <- run_dss(targets, args, cmp=cmp[[1]], smoothing=TRUE)

