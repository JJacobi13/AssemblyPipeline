library(seqinr)

args <- commandArgs(TRUE)

faFile <- read.fasta(args[1])
sortedLen <- sort(unlist(lapply(faFile,length)), decreasing=T)
summedLen <- cumsum(sortedLen)
png(args[2])
plot(summedLen, main="A 50 plot of the assembly", ylab="Cumulative sum of contig lengths", xlab="Contig index")
dev.off()