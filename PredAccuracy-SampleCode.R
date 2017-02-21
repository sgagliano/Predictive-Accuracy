#Sample R Code
#MyData.txt is a text file (either space or tab-delimited) with a header. It contains a list of genetic variants (one variant per line) with at least the following two labelled columns: cls (a 0/1 binary indicator: 1=hit and 0=non-hit) and score (contains the prediction value).

#Recceiver operator characteristic curve
pdf("ROC.pdf")
x<-read.table("MyData.txt", h=T, as.is=T)
library(ROCR)
pred<-prediction(x$score, x$cls) 
perf <- performance( pred, "tpr", "fpr" )
plot(perf, lwd=5)
abline(0,1,lty=3)
dev.off()
#display area under the curve
performance(pred, "auc")
#display positive predictive values
performance(pred, "ppv")
#display negative predictive values
performance(pred, "npv")

#Histogram
require(plotrix)
x<-read.table("MyData.txt", h=T, as.is=T)
hits<-subset(x, x$cls==1)
nonhits<-subset(x, x$cls==0)
l<-list(hits$score, nonhits$score)
#adjust the start and end position and bin increments below
bins<-seq(0,1, by=0.05)
pdf("Histogram.pdf")
multhist(l, freq=F, xlab="Predicted Value", breaks=bins, col=c("black","grey"))
legend("top", title="Classifier", c("Hits", "Non-hits"), pch=c(15, 15), col=c("black","grey"))
dev.off()

#Box plot
pdf("Boxplot.pdf")
x<-read.table("MyData.txt", h=T, as.is=T)
hits<-subset(x, x$cls=="1")
nonhits<- subset(x, x$cls=="0")
boxplot(hits$score, nonhits$score, xlab="Classification", ylab="Prediction", names=c("Hit","Non-hit"), ylim=c(0,1))
dev.off()

#Violin plot
library(vioplot) 
pdf("Violinplot.pdf")
x<-read.table("MyData.txt", h=T, as.is=T)
hits<-subset(x, x$cls=="1")
nonhits<- subset(x, x$cls=="0")
vioplot(hits$score, nonhits$score, names=c("Hit","Non-hit"), col="white", ylim=c(0,1))
title(xlab="Classification", ylab="Prediction")
dev.off()

#Quantile-quantile plot
pdf("Qqplot.pdf")
x<-read.table("MyData.txt", h=T, as.is=T)
hits<-subset(x, x$cls=="1")
nonhits<- subset(x, x$cls=="0")
qqplot(nonhits$score, hits$score, ylab="Hits", xlab="Non-hits", ylim=c(0,1), xlim=c(0,1))
abline(0,1, col="grey")
dev.off()

#Hypergeometric test
x<-read.table("MyData.txt", h=T, as.is=T)
hits<-subset(x, x$cls==1)
nonhits<-subset(x, x$cls==0)
res<-matrix(nrow=3,ncol=13)
row=1
col=0
BD<-length(nonhits[,1])
j<-length(hits[,1])
#prediction value bins ranging from less than 0.35 to between 0.9 and 0.95, increasing by increments of 0.5
for (i in seq(0.3,0.9,0.05))
{
col<-col+1
c<-length(subset(nonhits$score,nonhits$score<i+0.05 & nonhits$score>i))
a<-length(subset(hits$score, hits$score<i+0.05 & hits$score>i))
res[row,col]<-c/dim(nonhits)[1]
res[row+1,col]<-a/dim(hits)[1]
res[row+2,col]<-sum(phyper(a,j,BD-j,c, lower.tail=F))
}
#write a table to read in Excel
head<-c("p<0.35", "0.35<p<0.4", "0.4<p<0.45", "0.45<p<0.5", "0.5<p<0.55", "0.55<p<0.6","0.6<p<0.65", "0.65<p<0.7", "0.7<p<0.75","0.75<p<0.8", "0.8<p<0.85", "0.85<p<0.9","0.9<p<0.95")
table<-rbind(head, res)
write.table(table, "Hypergeometric.csv", sep=",", row.names=F, col.names=F, quote=F)
#the first row is the frequency of non-hits
#the second row is the frequency of the hits
#the third row is the hypergemoetric p-value

#Mann-Whitney U test
x<-read.table("MyData.txt", h=T, as.is=T)
nonhits<-subset(x, x$cls==0)
hits<-subset(x, x$cls==1)
wilcox.test(nonhits$score, hits$score)

#Asymptotic Generalized Cochran-Mantel-Haenszel Test
library("coin")
x<-read.table("MyData.txt", h=T, as.is=T)
nonhits<-subset(x, x$cls==0)
hits<-subset(x, x[$cls==1)
counts<-matrix(nrow=2,ncol=13)
row=1
col=0
for (i in seq(0.3,0.9,0.05))
{
col<-col+1
c<-length(subset(nonhits$score,nonhits$score<i+0.05 & nonhits$score>i))
a<-length(subset(hits$score, hits$score<i+0.05 & hits$score>i))
counts[row,col]<-c
counts[row+1,col]<-a
}
counts<-as.table(counts)
cmh_test(counts)
