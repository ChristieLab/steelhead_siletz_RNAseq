#=============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com
# Script created in version R 4.2.2 
# This script:
# Usage notes: 
#============================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
#setwd("~/files/", sep = "")) 


list.files()

#source("~/Drive/scripts.R", sep = ''))

dat <- read.table("dat.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
head(dat)

#=============================================================================================================#

main <- 386 # number of genes between main crosses
rec  <- 5   # number of genes between reciprocal crosses

total <- main + rec

OUT <- NULL
for(n in 1:1000000){
  sam1 <- sample(c(1,2), total, replace = TRUE)
  out  <- length(which(sam1 == 1))
  OUT  <- rbind(OUT, out)
  
}

range(OUT)

hist(OUT, xlab = "Number of DEGs", xlim = c(140, 390), ylim=c(0, 80000), col="plum2", breaks = 35, border="plum3", main="")
#abline(h=0)
abline(v=386, col="blue", lwd = 2, lty=2)



barplot(c(main, rec), ylim = c (0, 400), col="dodgerblue", cex.axis=1.25, width=2, space=1, xlim=c(1,8))

abline(h=0, lwd=2)




