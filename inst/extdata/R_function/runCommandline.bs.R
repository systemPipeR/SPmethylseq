##############################################################################
## Function to run NGS aligners including sorting and indexing of BAM files ##
##############################################################################
runCommandline_bs <- function(args, runid="01") {
	if(any(nchar(gsub(" {1,}", "", modules(args))) > 0)) {
	# if(system("module -V", ignore.stderr=TRUE)==1) { # Returns 1 if module system is present. This is a better solution, but run some test before committing it!
        for(j in modules(args)) moduleload(j) # loads specified software from module system
	}	
	commands <- sysargs(args)
	for(i in seq(along=commands)) {
		s <- system(commands[i])
        }
	
	return(s)
}

## Usage: 
# runCommandline(args=args)
