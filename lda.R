#!/usr/bin/Rscript --vanilla

source("http://aoki2.si.gunma-u.ac.jp/R/src/disc.R", encoding="euc-jp")

size = 100000

prof = read.csv("/home/hacone/src/AgIn/Profile/P6_profile_shuffled.dat", sep=",", header=F)
#data = prof[1:size,2:22]
#group = factor(prof[1:size,23])
data = prof[,2:22]
group = factor(prof[,23])
discans = disc(data, group)

print("LDA done")

svg("/home/hacone/LDA_P6.svg")
barplot(discans$d.function[1:21,1], names.arg=-10:10,
        xlab="relative position to CpG", ylab="coefficients of the LD function")
dev.off()

print("LDA.svg written")

ldavec = discans$d.function[1:21,1]

write(discans$d.function[1:21,1], file="/home/hacone/LDA_P6.dat", sep="\n")

print("all done")
