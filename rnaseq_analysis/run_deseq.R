#!/usr/bin/Rscript

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 4) {
      args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
      cat("
                      Arguments:
                      input - tab-delimitd file with counts
                      outputfilename - csv output file with results
		      outputfile2 - pdf with plots
		      fittype - local or parametric
		      N - for each sample 0 if untreated, 1 if treated  \n\n")
                                   
                                    q(save="no", status=1)
}
 
datafile <- args[1]
outputfile <- args[2]
outputfile2 <- args[3]
fittype <- args[4]

if ((fittype != "local") && (fittype != "parametric")) {
    cat(" Fit type can be only local or parametric ")
    q(save="no", status=1)
}

args <- args[5:length(args)]
countTable = read.table( datafile, sep="\t", header=TRUE, row.names=1 )
#print(countTable)

if ((length(args) != length(colnames(countTable))) || (sum(args == "0") + sum(args == "1") != length(args))) {
    cat(" You need to specify if sample is treated (1) or untreated(0) for every sample ")
    q(save="no", status=1)
}

f <- function(x) {
    f <- ifelse(x=="0", "untreated", "treated")
}
condition <- sapply(args, f)
n <- length(condition)


#load library, silently
suppressWarnings(suppressMessages(library( "DESeq" )))
title=paste("DESeq analysis of ", sum(condition=="untreated"), " untreated and ", sum(condition=="treated"), " treated samples."  ) 
print(Sys.time())
print(title)


pdf(outputfile2, paper='a4', title=title)
par(mfrow=c(2,1))
#DESeq
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )


if ( length(condition) > 3) {
    cds = estimateDispersions( cds, fitType=fittype )
} else {
    cds = estimateDispersions( cds, method="blind", sharingMode="fit-only", fitType="local" )
}
#print( fitInfo( cds ) )
par(oma=c(0,0,3,0), mar=c(5,5,3,5))
plotDispEsts( cds, main="Empirical (black dots) and fitted (red lines) dispersion values \nplotted against the mean of the normalised counts. " )

res = nbinomTest( cds, "untreated", "treated" )
plotMA(res, main="Plot of normalised mean versus log2 fold change \nfor the contrast untreated versus treated")
mtext(title, outer=TRUE, cex=1.5, line=1)

hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="Histogram of p-values")

par(xaxs='i', yaxs='i' )
plot.new()
s <- lapply(sizeFactors(cds), 
	   function(x) {
	    formatC(x, digits=4, format="g", flag="-") })
par(mar=c(0,0,0,0))
places=(0:(n-1))/n
text(0.1, 0.35, "Size Factors:", offset=0, pos=4)
text(places, 0.2, s, offset=0, pos=4)
text(0.1, 0.95, "Samples:", offse=0, pos=4)
text(places, 0.8, colnames(countTable), offset=0, pos=4) 
text(places, 0.65, condition, offset=0, pos=4) 
#text(0.1, 0.6, "Fit info:", offset=0, pos=4)
#text(0.2, 0.9, fData(cds), pos=4, offset=0)

#f <- fitInfo(cds)
#text(c(0.1), c(0.4), c(f$dispFunc), pos=4, offset=0)

garbage <- dev.off()

#outputfile<-file(outputfile)

res <- res[ order(res$pval), ] 
write.table(res, outputfile, append=FALSE, quote=FALSE, row.names = FALSE, sep="\t")
print(Sys.time())

q(save="no", status=0)
