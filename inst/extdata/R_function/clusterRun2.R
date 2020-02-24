###########################################################################################
## batchtools-based function to submit runCommandline jobs to queuing system of a cluster ##
###########################################################################################
## The advantage of this function is that it should work with most queuing/scheduling systems such as SLURM, Troque, SGE, ...
clusterRun2 <- function (args, FUN = runCommandline, conf.file = ".batchtools.conf.R",
                        template = "batchtools.slurm.tmpl", runid = "01", resourceList) {                                                                                                     
	## Validity checks of inputs
    if(class(args)!="SYSargs") stop("Argument 'args' needs to be assigned an object of class 'SYSargs'")
	if(class(FUN)!="function") stop("Value assigned to 'FUN' argument is not an object of class function.")
    if(!file.exists(conf.file)) stop("Need to point under 'conf.file' argument to proper config file. See more information here: https://mllg.github.io/batchtools/reference/makeRegistry.html.
                                     Note: in this file *.tmpl needs to point to a valid template file.")
    if(!file.exists(template)) stop("Need to point under 'template' argument to proper template file. Sample template files for different schedulers are available here: https://github.com/mllg/batchtools/blob/master/inst/templates/")
    ## batchtools routines

    f <- function(i, args, ...) FUN(args=args[i], ...)
    logdir1 <- paste0(normalizePath(results(args)), "/submitargs", runid, "_BJdb_", paste(sample(0:9, 4), collapse = ""))                                                                               
    reg <- makeRegistry(file.dir = logdir1, conf.file = ".batchtools.conf.R", packages = "systemPipeR")                                                                                                 
    ids <- batchMap(fun = f, seq(along = args), more.args = list(args = args, runid = runid), reg=reg)                                                                                                  
    done <- submitJobs(ids=ids, reg=reg, resources = resourceList)                                                                                                                                      
    return(reg)                                                                                                                                                                                             
}
## Usage: 
# resources <- list(walltime=3600, nodes=paste0("1:ppn=", cores(args)), memory=1024)
# reg <- clusterRun(args, conf.file = ".batchtools.conf.R", template = "batchtools.slurm.tmpl, runid="01", resourceList=resources)
# waitForJobs(reg=reg)
# getStatus(reg=reg)  
