############################################################
##  Usage: plot_insertsizes.R <input file> <output file>  ##
##  Input file is a histogram of insert sizes             ##
##  Output file is a plot of this histogram               ##
############################################################
args <- commandArgs(TRUE)
#read the given histogram (first command-line argument) into a table
dat=read.table(args[1])
#Create values for drawing a smooth line, to remove outliers from the plot.
smth <- smooth.spline(dat[,2], spar=0.5)$y
#calculate the mean insert size for better boundaries of the plot
meanInsertSize <- dat[which.max( dat[,2] ),1]
#add the smooth line to the data as a new column
dat[,3] <- smth
#Plot the histogram
png(args[2])
plot(dat[,1],dat[,2],
     main="Insert size distribution", 
     ylab="Number of reads", 
     xlab="Insert size",
     xlim=c(min(dat[,1]),meanInsertSize + (meanInsertSize-min(dat[,1]))),
     ylim=c(0,max(smth)))
#draw a red line of the smoothed line
lines(dat[,1],dat[,3], col="red")
dev.off()