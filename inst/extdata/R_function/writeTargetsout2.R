##############################################
## Additional utilities for SYSargs objects specific to BS-seq ##
##############################################

## Convenience write function for targetsout(args)
writeTargetsout2 <- function(x, file="default", silent=FALSE, overwrite=FALSE, folder=FALSE, ...) {
  if(class(x)!="SYSargs") stop("x needs to be 'SYSargs' object")
  targets <- targetsout(x)
 
  ## If the targetsout(x) are a folder, it is necessary create a targets file with the bam files.
  if(folder==TRUE){
  bam.file <- as.list(targets$FileName)
  ## Collect the information
  list <- list()
  for (i in seq(along=bam.file)){
    file.list <- as.list(list.files(paste0(bam.file[i]), "*.bam", full.names=TRUE))
    list[length(list)+1] <- file.list
  }
  targets$FileName <- list
  }
  
  if(file=="default") {
    file <- paste("targets_", software(x), ".txt", sep="")
    file <- gsub(" {1,}", "_", file)
  } else {
    file <- file
  }
  if(file.exists(file) & overwrite==FALSE) stop(paste("I am not allowed to overwrite files; please delete existing file:", file, "or set 'overwrite=TRUE'"))
  headerlines <- targetsheader(x)	
  targetslines <- c(paste(colnames(targets), collapse="\t"), apply(targets, 1, paste, collapse="\t"))
  writeLines(c(headerlines, targetslines), file, ...)
  if(silent!=TRUE) cat("\t", "Written content of 'targetsout(x)' to file:", file, "\n")
}
## writeTargetsout(x=args, file="default") 