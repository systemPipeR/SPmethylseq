##############################################################################
## Function to run NGS aligners including sorting and indexing of BAM files ##
##############################################################################
runCommandline2 <- function(args, runid="01", make_bam=FALSE, ...) {
	if(any(nchar(gsub(" {1,}", "", modules(args))) > 0)) {
	# if(system("module -V", ignore.stderr=TRUE)==1) { # Returns 1 if module system is present. This is a better solution, but run some test before committing it!
        for(j in modules(args)) moduleload(j) # loads specified software from module system
	}	
	commands <- sysargs(args)
	completed <- file.exists(outpaths(args))
	names(completed) <- outpaths(args)
	logdir <- results(args)
	for(i in seq(along=commands)) {
		## Run alignmets only for samples for which no BAM file is available.
		if(as.logical(completed)[i]) {
			next()
		} else {
			## Create soubmitargsID_command file
			cat(commands[i], file=paste(logdir, "submitargs", runid, sep=""), sep = "\n", append=TRUE)
        		## Run executable 
			command <- gsub(" .*", "", as.character(commands[i]))
			commandargs <- gsub("^.*? ", "",as.character(commands[i]))
			
			## Execute system command; note: BWA needs special treatment in stderr handling since it writes 
                        ## some stderr messages to sam file if used with system2()
			if(software(args) %in% c("bwa aln", "bwa mem")) {
				stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=FALSE)
			} else if(software(args) %in% c("bash_commands", "samtools")) {
				stdout <- system(paste(command, commandargs))
			} else {
				stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=TRUE)
			}
			
			## Create submitargsID_stdout file
			cat(commands[i], file=paste(logdir, "submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
			cat(unlist(stdout), file=paste(logdir, "submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
			## Conditional postprocessing of results
			if(make_bam==TRUE) {
                if(grepl(".sam$", outfile1(args)[i])) { # If output is *.sam file (e.g. Bowtie2)
				    asBam(file=outfile1(args)[i], destination=gsub("\\.sam$", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)
				    unlink(outfile1(args)[i])
			    } else if(grepl("vcf$|bcf$|xls$|bed$", outpaths(args)[i])) {
                    dump <- "do nothing"
			    } else { # If output is unindexed *.bam file (e.g. Tophat2)
				    sortBam(file=names(completed[i]), destination=gsub("\\.bam$", "", names(completed[i])))
        			        indexBam(names(completed[i]))
			    }
		    }
        }
	}
	bamcompleted <- gsub("sam$", "bam$", file.exists(outpaths(args)))
	names(bamcompleted) <- SampleName(args)
	cat("Missing alignment results (bam files):", sum(!as.logical(bamcompleted)), "\n"); cat("Existing alignment results (bam files):", sum(as.logical(bamcompleted)), "\n")
	return(bamcompleted)
}

## Usage: 
# runCommandline(args=args)
