read.align.bs <- function(path, plot=TRUE, table=TRUE){
  ## Find the path
  if(missing(path)) { 
    path <- getwd()
    path <- paste(path, "/results",sep="")
  } else { 
    path <- path
  }
  ## Read the files
  file.list <- as.list(list.files(path, "*.out", full.names=TRUE))
  colnames <- as.list(list.files(paste0(file.list[1]), "*.txt", full.names=TRUE))
  colnames <- read.table(paste0(colnames[1]), header=F, row.names=NULL, sep="\t", fill=T, quote="\"") 
  colnames <- colnames[c(6,8,10,12,14,16,19,22,25,28,35,37,39,41,45,47,49,51,53,55,57,59), , drop=F]
  align_results <- data.frame(matrix(nrow = 22));rownames(align_results) <- colnames[,1:1] ##TODO: remove ":" rownames
  ## Collect the information
  for (i in seq(along=file.list)){
    file.list1 <- as.list(list.files(paste0(file.list[i]), "*.txt", full.names=TRUE))
    align <- read.table(paste0(file.list1[1]), header=F, row.names=NULL, sep="\t", fill=T, quote="\"")
    align_sub <- align[c(7,9,11,13,15,17,20,23,26,36,38,40,42,44,46,48,50,52,54,56,58,60), , drop=F]
    align_results <- cbind(align_results, align_sub)
  }
  align_results <- align_results[,-1]
  colnames <- gsub("^.*/([A-z0-9]*).fastq_.*","\\1", file.list)
  colnames(align_results) <- colnames
  ## Plot Methylation levels
  if(plot==TRUE){
  sub_a <- align_results[19:21,]
  names <- c("CpG", "CHG", "CHH")
  rownames(sub_a) <- names
  
  CpG <- data.frame(Samples=colnames(sub_a), Methylation=c(t(sub_a[1:1,])), Type=(rep("CpG")))
  CHG <- data.frame(Samples=colnames(sub_a), Methylation=c(t(sub_a[2:2,])), Type=(rep("CHG")))
  CHH <- data.frame(Samples=colnames(sub_a), Methylation=c(t(sub_a[3:3,])), Type=(rep("CHH")))
  df_methyl <- rbind(CHH, CHG, CpG)
  
  p_align <- ggplot(data = df_methyl, aes(x=Samples, y=Methylation, fill=Type)) +
    geom_col(aes(fill = Type)) + ylab("Methylation (%)") +
    geom_text(aes(label = Methylation), position = position_stack(vjust = 0.5), color="white", size=3.5) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "grey"))  
    print (p_align)
  ## Plot alignment results
  sub_b <- align_results[c(1:2,4:6),]
  names_b <- c("Total_Sequences", "Unique_Alignments1", "No_Alignments1", "Multiple_Alignments1", "No_Genomic_Sequence1") 
  sub_b <- data.frame(t(sub_b))
  colnames(sub_b) <- names_b
  Unique_Alignments <-  data.frame(Samples=rownames(sub_b), Percent=(as.numeric(as.character(sub_b$Unique_Alignments1))/as.numeric(as.character(sub_b$Total_Sequences)))*100, Type=(rep("Unique Alignments")))
  No_Alignments <- data.frame(Samples=rownames(sub_b), Percent=(as.numeric(as.character(sub_b$No_Alignments1))/as.numeric(as.character(sub_b$Total_Sequences)))*100, Type=(rep("No Alignments")))
  Multiple_Alignments <- data.frame(Samples=rownames(sub_b), Percent=(as.numeric(as.character(sub_b$Multiple_Alignments1))/as.numeric(as.character(sub_b$Total_Sequences)))*100, Type=(rep("Multiple Alignments")))
  No_Genomic_Sequence <- data.frame(Samples=rownames(sub_b), Percent=(round(as.numeric(as.character(sub_b$No_Genomic_Sequence1))/as.numeric(as.character(sub_b$Total_Sequences))))*100, Type=(rep("No Genomic Sequence")))
  
  df_align <- rbind(No_Alignments, Multiple_Alignments, Unique_Alignments)
  p_methylation <- ggplot(data = df_align, aes(x=Samples, y=Percent, fill=Type)) +
    geom_col(aes(fill = Type)) + ylab("Alignment (%)") +
    geom_text(aes(label = round(Percent, 2)), position = position_stack(vjust = 0.5), color="white", size=3.5)+
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "grey"))
  ## Save plots
  pdf("results/align.bs_plot.pdf")
  print (p_align)
  dev.off()
  pdf("results/methyl.levels_plot.pdf")
  print (p_methylation)
  dev.off()
}
  if(table==TRUE){
    write.table(align_results, "./results/align.bs_results.xls", quote=FALSE, sep="\t", row.names=TRUE)
  }
  return(align_results)
}

## Usage:
# read_align <- read.align.bs()
