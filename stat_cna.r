genQuantile <- function(dat, type, main=NULL){
  library(ggplot2)

  # type v: indicates that the data is directly from VarScan
  # type vc: indicates that the data is from VarScan + copyCaller
  # type vcc: indicates that the data is from VarScan + copyCaller + cbs
  title = ""
  if (type=='v'){
    title <- paste("The 10-quantiles of the Segment Lengths (VarScan)", sep="")
  }
  if (type=="vc"){
    title <- paste("The 10-quantiles of the Segment Lengths (VarScan + copyCaller)", sep="")
  }
  if (type=="vcc"){
    title <- paste("The 10-quantiles of the Segment Lengths in Whole-exome Space (VarScan + copyCaller + cbs", sep="")
  }
  if (is.null(main)==FALSE){
    title=main
  }
  print(title,justify=1)

  p <- c(.01,.03,.05,.1,.15,.20,.50,.70,.80,.90)
  print(quantile(dat, probs=p))


  # plot quantiles with ggplot2
  # http://stackoverflow.com/questions/8894617/plot-quantiles-with-ggplot2
  # http://stackoverflow.com/questions/2359723/how-to-add-a-title-to-a-ggplot-when-the-title-is-a-variable-name
  dat1 <- data.frame(q=quantile(dat, probs=p), prob=p)

  ggplot(aes(x=prob, y=q), data=dat1)+geom_line()+ggtitle(title)+xlab('Quantile')+ ylab('Segment Lengths')+ggtitle("10-Quantiles of Segment Length in VarScan+cbs Output in Sample 1767")

}

genHist <- function(dat, title){
  # Take the cna data (which is from read.table() directly)
  colnames(dat)[4] <- "num_positions" # I did this because the 4th column name of cbs output is num.mark. I need to change it to num_position to keep it uniform with other data frame.
  ggplot(dat, aes(x=num_positions))+geom_histogram(binwidth=10) + scale_y_log10() + ggtitle(title) + xlab("Segment Length") + ylab("Log Count")
  #+ylim(c(0,max(dat$num_positions)))
}

genBoxplot <- function(dat){
  ggplot(dat, aes(x=0,y=num_positions))+geom_boxplot()+scale_y_log10()+ggtitle("Boxplot of log lengths")
}


plotAll <- function(files=c(1767)){

  # files=c(1767, 1952, 1954, 2036, 2241, 2277, 2295, 2405, 2668, 2825)
  bases <- files
  print(bases)
  bases <- as.character(bases)

  # v is the filenames of varscan output.
  v <- paste(bases, ".copynumber", sep="")
  print(paste("I will take the following files as VarScan copynumber output", paste(v, collapse=", ")))
  # vc is the filenames of varscan+copyCaller output
  vc <- paste(v, ".called", sep="")
  print(paste("I will take the following files as VarScan copynumber + copyCaller output", paste(vc, collapse=", ")))
  # vcc is the filenames of varscan+copyCaller+cbs output
  vcc <- paste(vc, ".cbs", sep="")
  print(paste("I will regard the following files as the output of VarScan copynumber + copyCaller + cbs", paste(vcc, collapse=", ")))

  # Split the page into several sections to plot different graphs

  dfs = list()
  len <- length(bases)
  # Read all the files into data frames
  for (i in 1:len){
    dfs[[i]] = read.table(v[i], header=TRUE)
  }

##################################
  # Quantile
  for (i in 1:len){
    genQuantile(dfs[[i]]$num_positions,'v')
  }
  
  # Hist (log)
  for (i in 1:len){
    dev.new()
    title=paste("Distribution of segment length in ", v[i], sep="")
    genHist(dfs[[i]], title)
  }

  
# Since the boxplot in ggplot requires data from ONE dataframe, the 
# following code will concatenate all the cbs file into one file with
# their sample name
  # Boxplot
  cbs1767.df <- read.table("1767.copynumber.called.cbs", header=TRUE)
  cbs1952.df <- read.table("1952.copynumber.called.cbs", header=TRUE)
  cbs1954.df <- read.table("1954.copynumber.called.cbs", header=TRUE)
  cbs2036.df <- read.table("2036.copynumber.called.cbs", header=TRUE)
  cbs2241.df <- read.table("2241.copynumber.called.cbs", header=TRUE)
  cbs2277.df <- read.table("2277.copynumber.called.cbs", header=TRUE)
  cbs2295.df <- read.table("2295.copynumber.called.cbs", header=TRUE)
  cbs2405.df <- read.table("2405.copynumber.called.cbs", header=TRUE)
  cbs2668.df <- read.table("2668.copynumber.called.cbs", header=TRUE)
  cbs2825.df <- read.table("2825.copynumber.called.cbs", header=TRUE)

  cbs1767.df[["name"]] = rep("1767", times=length(cbs1767.df[[1]]))
  cbs1952.df[["name"]] = rep("1952", times=length(cbs1952.df[[1]]))
  cbs1954.df[["name"]] = rep("1954", times=length(cbs1954.df[[1]]))
  cbs2036.df[["name"]] = rep("2036", times=length(cbs2036.df[[1]]))
  cbs2241.df[["name"]] = rep("2241", times=length(cbs2241.df[[1]]))
  cbs2277.df[["name"]] = rep("2277", times=length(cbs2277.df[[1]]))
  cbs2295.df[["name"]] = rep("2295", times=length(cbs2295.df[[1]]))
  cbs2405.df[["name"]] = rep("2405", times=length(cbs2405.df[[1]]))
  cbs2668.df[["name"]] = rep("2668", times=length(cbs2668.df[[1]]))
  cbs2825.df[["name"]] = rep("2825", times=length(cbs2825.df[[1]]))

  write.table(cbs1767.df, file="1767.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs1952.df, file="1952.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs1954.df, file="1954.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs2036.df, file="2036.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs2241.df, file="2241.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs2277.df, file="2277.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs2295.df, file="2295.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs2405.df, file="2405.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs2668.df, file="2668.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")
  write.table(cbs2825.df, file="2825.cbs.group", row.name=FALSE, col.name=FALSE, quote=F, sep="\t")

  all.dat <- read.table("all.cbs.group", header=FALSE)
  ggplot(all_dat, aes(x=name, y=num.mark))+geom_boxplot()
  ggplot(all.dat, aes(x=V6, y=V4)) + geom_boxplot() + scale_y_log10() + xlab("sample") + ylab("log segment lengths in whole-exome space")+ggtitle("Segment Lengths by Sample")

  all.dat[['V7']] = all.dat[['V3']] - all.dat[['V2']] # On genome-scale
  ggplot(all.dat, aes(x=V6, y=V7)) + geom_boxplot() + scale_y_log10() + xlab("sample") + ylab("log segment lengths in whole-genome space")+ggtitle("Segment Lengths by Sample")
}
