library(seqinr)

args <- commandArgs(TRUE)

faFile <- read.fasta(args[1])
sortedLen <- sort(unlist(lapply(faFile,length)), decreasing=T)
summedLen <- cumsum(sortedLen)
summedLen <- summedLen/1000000
png(args[2])
plot(summedLen, type="l", main="A 50 plot of the assembly", ylab="Cumulative sum of sequence lengths (MB)", xlab="Sequence index")
dev.off()