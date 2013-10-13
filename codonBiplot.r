#!/usr/bin/env Rscript

progname = commandArgs(FALSE)[4]
args = commandArgs(TRUE)

unchangedMatrix = args[1]
novelMatrix = args[2]

um = read.csv(unchangedMatrix, sep="\t", header=F)
nm = read.csv(novelMatrix, sep="\t", header=F)

pdf("CodonUsage.pca.pdf",width=8.5,height=11)
upca = prcomp(um,scale=TRUE)
biplot(upca,var.axes=F,xlabs=rep('*',nrow(um)),ylabs=rep('.',ncol(um)), col=c('black','white'))
p = predict(upca,nm)
points(p[,1],p[,2],col='blue',pch=16)
