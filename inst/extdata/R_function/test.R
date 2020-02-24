library(methylKit) 

file.list2=list(system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"), system.file("extdata", "test.fastq_bismark.sorted.min.sam", 
               package = "methylKit"), system.file("extdata", "test.fastq_bismark.sorted.min.sam", 
              package = "methylKit"), system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                package = "methylKit"))



## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(dpi = 75)
knitr::opts_chunk$set(cache = FALSE)

## ----eval=TRUE,echo=FALSE,warning=FALSE,message=FALSE--------------------
#devtools::load_all(".")

## ---- eval=TRUE, echo=FALSE----------------------------------------------
tab <- read.table( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),header=TRUE,nrows=5)
tab
#knitr::kable(tab)

## ----message=FALSE-------------------------------------------------------
library(methylKit)
file.list=list( system.file("extdata", 
                            "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata",
                            "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", 
                            "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", 
                            "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),
           assembly="hg18",
           treatment=c(1,1,0,0),
           context="CpG"
           )

myobj2 <- methRead(a,
           sample.id=list("C1","C2","C3","S1", "S2", "S3"),
           assembly="tair10",
           treatment=c(0,0,0,1,1,1),
           context="CpG"
           )


#
# ---- message=FALSE,warning=FALSE----------------------------------------
library(methylKit)
file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawListDB object: myobjDB 
# and save in databases in folder methylDB
myobjDB=methRead(file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),
           assembly="hg18",
           treatment=c(1,1,0,0),
           context="CpG",
           dbtype = "tabix",
           dbdir = "methylDB"
           )

print(myobjDB[[1]]@dbpath)



## ---- eval=FALSE---------------------------------------------------------
#  my.methRaw=processBismarkAln( location =
#                                  system.file("extdata",
#                                                  "test.fastq_bismark.sorted.min.sam",
#  	                                              package = "methylKit"),
#                           sample.id="test1", assembly="hg18",
#                           read.context="CpG", save.folder=getwd())

## ------------------------------------------------------------------------
getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)

## ------------------------------------------------------------------------
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

## ------------------------------------------------------------------------
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

## ------------------------------------------------------------------------
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

## ------------------------------------------------------------------------
meth=unite(myobj, destrand=FALSE)

## ------------------------------------------------------------------------
head(meth)

## ----eval=FALSE----------------------------------------------------------
#  # creates a methylBase object,
#  # where only CpGs covered with at least 1 sample per group will be returned
#  
#  # there were two groups defined by the treatment vector,
#  # given during the creation of myobj: treatment=c(1,1,0,0)
#  meth.min=unite(myobj,min.per.group=1L)

## ------------------------------------------------------------------------
getCorrelation(meth,plot=TRUE)

## ------------------------------------------------------------------------
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

## ----message=FALSE-------------------------------------------------------
hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)

## ------------------------------------------------------------------------
PCASamples(meth, screeplot=TRUE)

## ------------------------------------------------------------------------
PCASamples(meth)

## ------------------------------------------------------------------------
# make some batch data frame
# this is a bogus data frame
# we don't have batch information
# for the example data
sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),
                            age=c(19,34,23,40))

as=assocComp(mBase=meth,sampleAnnotation)
as

# construct a new object by removing the first pricipal component
# from percent methylation value matrix
newObj=removeComp(meth,comp=1)

## ------------------------------------------------------------------------
mat=percMethylation(meth)

# do some changes in the matrix
# this is just a toy example
# ideally you want to correct the matrix
# for batch effects
mat[mat==100]=80
 
# reconstruct the methylBase from the corrected matrix
newobj=reconstruct(mat,meth)

## ----warning=FALSE-------------------------------------------------------
tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000)
head(tiles[[1]],3)

## ------------------------------------------------------------------------
myDiff=calculateDiffMeth(meth)

## ------------------------------------------------------------------------
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#
#
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

## ------------------------------------------------------------------------
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

## ----eval=FALSE----------------------------------------------------------
#  
#  sim.methylBase1<-dataSim(replicates=6,sites=1000,
#                           treatment=c(rep(1,3),rep(0,3)),
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
#                          )
#  
#  my.diffMeth<-calculateDiffMeth(sim.methylBase1[1:,],
#                                  overdispersion="MN",test="Chisq",mc.cores=1)

## ----eval=FALSE----------------------------------------------------------
#  
#  covariates=data.frame(age=c(30,80,34,30,80,40))
#  sim.methylBase<-dataSim(replicates=6,sites=1000,
#                          treatment=c(rep(1,3),rep(0,3)),
#                          covariates=covariates,
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
#                          )
#  my.diffMeth3<-calculateDiffMeth(sim.methylBase,
#                                  covariates=covariates,
#                                  overdispersion="MN",test="Chisq",mc.cores=1)

## ---- eval=FALSE---------------------------------------------------------
#  myDiff=calculateDiffMeth(meth,mc.cores=2)

## ------------------------------------------------------------------------
library(genomation)

# read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                           package = "methylKit"))
#
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
#
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

## ------------------------------------------------------------------------
# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                        package = "methylKit"),
                           feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")

## ------------------------------------------------------------------------
promoters=regionCounts(myobj,gene.obj$promoters)

head(promoters[[1]])

## ------------------------------------------------------------------------
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

# target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))

## ------------------------------------------------------------------------
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

## ------------------------------------------------------------------------
plotTargetAnnotation(diffAnn,precedence=TRUE,
    main="differential methylation annotation")

## ------------------------------------------------------------------------
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
       main="differential methylation annotation")

## ------------------------------------------------------------------------
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

## ------------------------------------------------------------------------
class(meth)
as(meth,"GRanges")
class(myDiff)
as(myDiff,"GRanges")

## ------------------------------------------------------------------------
class(myobjDB[[1]])

## ----eval=FALSE----------------------------------------------------------
#  as(myobjDB[[1]],"methylRaw")

## ----eval=FALSE----------------------------------------------------------
#  data(methylKit)
#  
#  objDB=makeMethylDB(methylBase.obj,"exMethylDB")
#  

## ------------------------------------------------------------------------
select(meth,1:5) # get first 10 rows of a methylBase object
myDiff[21:25,] # get 5 rows of a methylDiff object

## ----message=FALSE,warning=FALSE,eval=FALSE------------------------------
#  library(GenomicRanges)
#  my.win=GRanges(seqnames="chr21",
#                 ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )
#  
#  # selects the records that lie inside the regions
#  selectByOverlap(myobj[[1]],my.win)

## ----eval=FALSE----------------------------------------------------------
#  # creates a new methylRawList object
#  myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
#  # creates a new methylBase object
#  meth2 =reorganize(meth,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )

## ---- eval=FALSE---------------------------------------------------------
#  # creates a matrix containing percent methylation values
#  perc.meth=percMethylation(meth)

## ---- eval=FALSE---------------------------------------------------------
#   download.file("https://dl.dropboxusercontent.com/u/1373164/H1.chr21.chr22.rds",
#                 destfile="H1.chr21.chr22.rds",method="curl")
#  
#   mbw=readRDS("H1.chr21.chr22.rds")
#  
#   # it finds the optimal number of componets as 6
#   res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10)
#  
#   # however the BIC stabilizes after 4, we can also try 4 componets
#   res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:4)
#  
#   # get segments to BED file
#   methSeg2bed(res,filename="H1.chr21.chr22.trial.seg.bed")
#  
#  

## ------------------------------------------------------------------------
sessionInfo() 

## ----eval=TRUE,echo=FALSE------------------------------------------------
# tidy up                  
rm(myobjDB)              
unlink(list.files(pattern = "methylDB",full.names = TRUE),recursive = TRUE)

