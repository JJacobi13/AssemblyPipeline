######################################################
## Makes Kmer plot from jellyfish histogram         ##
## First argument = histogram from jellyfish        ##
## Second argument = output file                    ##
######################################################
args <- commandArgs(TRUE) 
#Read file into dataframe
dat <- read.table(args[1], row.names=NULL)
#The number of kmers is the number of occurrences * counts of kmers with this occurence
dat <- dat[,1] * dat[,2]
#Remove first line (and empty lines before)
dat <- dat[3:length(dat)-1]
#Create a smooth path through all datapoints to retrieve only a single peak
smth <- smooth.spline(dat, spar=0.5)$y
#Calculate the peak of the smooth line
peak <- which(diff(sign(diff(smth)))==-2)+1
print(peak)

#Draw the plot and write to the output file
png(args[2])
if(length(peak) == 1){
  plot(dat, main="kmer coverage histogram",ylab="Volume of kmers", xlab="Counts of a k-mer",ylim=c(0,dat[peak]*1.1),xlim=c(2,peak*2))
}else{
  plot(dat, main="kmer coverage histogram",ylab="Volume of kmers", xlab="Counts of a k-mer",ylim=c(0,dat[peak[1]]*1.1))
}
#lines(smth, col="blue") #Show the smooth line which determines the peak
#Draw an arrow to the calculated peak
arrows(peak,0, peak,dat[peak]*0.98)
dev.off()
