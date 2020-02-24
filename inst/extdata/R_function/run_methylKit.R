####################################################################
## Import Bismark results with methylKit package ##
####################################################################
import.bs.data <- function(targets, args, cmp=cmp[[1]], assembly="TAIR10"){
  ##Settings
  if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
  samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
  targets_list <- as.list(targetsout(args)$FileName)
  report <- list()
  for (i in seq(along=targets_list)){
    file.list <- as.list(list.files(paste0(targets_list[i]), "*.cov", full.names=TRUE))
    report[length(report)+1] <- file.list
    }
  names(report) <- names(samples)
  meth <- list()
  ## Loop for each comparisons of the analysis
  for (i in seq(along=cmp[,1])) {
    treat0 <- rep(0, length(samples[samples %in% cmp[i,][1]]))
    treat1 <- rep(1, length(samples[samples %in% cmp[i,][2]]))
    sub <- c(report[names(report) %in% names(samples[samples %in% cmp[i,][1]])], report[names(report) %in% names(samples[samples %in% cmp[i,][2]])])
    ## import data
    methRaw <- methylKit::methRead(location=sub, sample.id=as.list(names(sub)), assembly=assembly, 
                                   dbtype = NA, pipeline = "bismarkCoverage", header = FALSE, skip = 0, sep = "\t",
                                   context = "CpG", resolution = "base", treatment=c(treat0, treat1), dbdir = getwd(), mincov = 10)
    meth <- c(meth, list(methRaw))
    }
  names <- apply(cmp, 1, paste, collapse="-")
  names(meth) <- names
  return(meth)
}

## Usage:
# targetspath <- "targets.txt"
# targets <- read.delim(targetspath, comment="#")
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# args <- systemArgs(sysma="param/bismark_methyl_extractor.param",
#                    mytargets="targets_bismark.txt")
# import.bs <- import.bs.data(targets, args, cmp=cmp[[1]])


####################################################################
## Run methylKit with cytosine methylation report from Bismark ##
####################################################################
run_methylKit <- function(targets, args, cmp=cmp[[1]], assembly="TAIR10", filter, plot=TRUE){
  ##Settings
  if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
  samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
  targets_list <- as.list(targetsout(args)$FileName)
  report <- list()
  for (i in seq(along=targets_list)){
    file.list <- as.list(list.files(paste0(targets_list[i]), "*.cov", full.names=TRUE))
    report[length(report)+1] <- file.list
  }
  names(report) <- names(samples)
  resultslist <- list()
  for (i in seq(along=cmp[,1])) {
    treat0 <- rep(0, length(samples[samples %in% cmp[i,][1]]))
    treat1 <- rep(1, length(samples[samples %in% cmp[i,][2]]))
    sub <- c(report[names(report) %in% names(samples[samples %in% cmp[i,][1]])], report[names(report) %in% names(samples[samples %in% cmp[i,][2]])])
    ## import data ##TODO: Future, if the user provides the methRaw data, we can skip this step
    methRaw <- methylKit::methRead(location=sub, sample.id=as.list(names(sub)), assembly=assembly, 
                                   dbtype = NA, pipeline = "bismarkCoverage", header = FALSE, skip = 0, sep = "\t",
                                   context = "CpG", resolution = "base", treatment=c(treat0, treat1), dbdir = getwd(), mincov = 10)
    
    ## Merging samples between two groups
    merged_samples <- unite(methRaw, destrand = FALSE)
    write.table(merged_samples, paste0("./results/merge_samples_", paste0(cmp[i,][1], sep="-", cmp[i,][2]), ".xls"), row.names=FALSE, quote=FALSE, sep="\t")
    print(paste0("Written the file: ", paste0("./results/merge_samples_", paste0(cmp[i,][1], sep="-", cmp[i,][2]), ".xls")))
    ## Finding differentially methylated bases or regions
    suppressWarnings(dmr_M <- calculateDiffMeth(merged_samples))
    # get differentially methylated bases
    dmr_M_filter <- getMethylDiff(dmr_M, qvalue=filter["qvalue"], difference=filter["difference"])
    resultslist <- c(resultslist, list(dmr_M_filter))
    ## Plot
    if(plot==TRUE){
    pdf(paste0("./results/DMR_methylkit_", paste0(cmp[i,][1], sep="-", cmp[i,][2]), ".pdf"))
    diffMethPerChr(dmr_M, plot=T, qvalue.cutoff=filter["qvalue"], meth.cutoff=filter["difference"])
    dev.off()
    }
  }
  names <- apply(cmp, 1, paste, collapse="-")
  names(resultslist) <- names
  return(resultslist)
}


## Usage:
# targetspath <- "targets.txt"
# targets <- read.delim(targetspath, comment="#")
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# args <- systemArgs(sysma="param/bismark_methyl_extractor.param",
#                    mytargets="targets_bismark.txt")
# dmr.DF <- run_methylKit(targets, args, cmp=cmp[[1]], assembly="TAIR10", filter=c(qvalue=0.01, difference=25), plot=TRUE)



