.libPaths("C:/Users/LCOIN/R-4.0.2/library")

library('jsonlite');
library('ggplot2')
library(tidyr)
library(reshape2)
RHOME="../../R"
#setwd("../data/new_samples")
source(paste(RHOME, "plotPlasmaFuncs.R", sep="/"));

names = c("1524", "1249", "1494", "1084","065", "098");




#x = 2:1999
x = 2:1999


x = 2:1999
window=1
bins = NULL
P0 = readDataAll(names[1:4], '_P0_read_length_count.json', x,"Tumour_all",window=window, bins = bins)
shared = readDataAll(names[1:4], '_P0.shared.somatic_reads_read_length_count.json', x,"Tumour_shared",window=window, bins = bins) #'_P0_shared_read_length_count.json', x)
unique = readDataAll(names[1:4], '_P0.unique.somatic_reads_read_length_count.json', x, "Tumour_unique",window=window, bins = bins)
P01 = readDataAll(names[5:6], '_P0_read_length_count.json', x, "Benign_all",window=window, bins = bins)
merged = readDataAll(names[5:6], '_P0.somatic_reads_read_length_count.json', x, "Benign_filtered",window=window, bins = bins)

mergedT = unique + shared
#mergedT = readDataAll(names[1:4], '_P0.somatic_reads_read_length_count.json', x,"Tumour_filtered",window=window, bins = bins)
x = as.numeric(dimnames(P0)[[1]])

.cumsum1<-function(x) {
  cum = cumsum(x)
  tot = sum(x)
  tot -cum
}

P0_c = data.frame(cbind(apply(P0,2,cumsum), apply(P0,2,.cumsum1)))
P01_c = data.frame(cbind(apply(P01,2,cumsum), apply(P01,2,.cumsum1)))
merged_c = data.frame(cbind(apply(merged,2,cumsum), apply(merged,2,.cumsum1)))

shared_c = data.frame(cbind(apply(shared,2,cumsum), apply(shared,2,.cumsum1)))
unique_c = data.frame(cbind(apply(unique,2,cumsum), apply(unique,2,.cumsum1)))
mergedT_c =  data.frame(cbind(apply(mergedT,2,cumsum), apply(mergedT,2,.cumsum1)))

unique_c[160,]+shared_c[160,] - mergedT_c[160,] ## should be all 0


ncol = dim(P0_c)[2]
ncol1 = dim(P01_c)[2]
ratios = matrix(nrow = dim(P0_c)[1],ncol = ncol)
ratios_1 = matrix(nrow = dim(P0_c)[1],ncol = ncol)
ratios_2 =  matrix(nrow = dim(P0_c)[1],ncol = ncol)
ratios_4 =  matrix(nrow = dim(P0_c)[1],ncol = ncol)
ratios_5 =  matrix(nrow = dim(P0_c)[1],ncol = ncol)
ratios_3 =  matrix(nrow = dim(P01_c)[1],ncol = ncol1)
for(i in 1:ncol1){
  ratios_3[,i] = merged_c[,i]/P01_c[,i]
}
for(i in 1:ncol){
  ratios[,i] = shared_c[,i]/mergedT_c[,i]
  ratios_1[,i] = mergedT_c[,i]/P0_c[,i]
  ratios_2[,i] = unique_c[,i]/mergedT_c[,i]
  ratios_4[,i] = shared_c[,i]/P0_c[,i]
  ratios_5[,i] = shared_c[,i]/unique_c[,i]
  
}
logx=T
xlim = c(0,1000)
allr = list("Shared to filtered"=ratios, "Tumour filtered to total" = ratios_1,
            "Unique to filtered"=ratios_2,
            "Shared to total"=ratios_4, "Shared to unique" = ratios_5)
allr1 = list( "Benign filtered to total" = ratios_3)
for(i in 1:length(allr)){
ggps =.convertRatios(x,allr[[i]],names(P0),main=names(allr)[[i]],logx=logx, xlim = xlim)
ggsave(paste(names(allr)[i],1,".png",sep=""), plot=ggps[[1]], width = 30, height = 30, units = "cm")
ggsave(paste(names(allr)[i],2,".png",sep=""), plot=ggps[[2]], width = 30, height = 30, units = "cm")
}
for(i in 1:length(allr1)){
  ggps =.convertRatios(x,allr1[[i]],names(P01),main=names(allr1)[[i]],logx=logx, xlim = xlim)
  ggsave(paste(names(allr1)[i],1,".png",sep=""), plot=ggps[[1]], width = 30, height = 30, units = "cm")
  ggsave(paste(names(allr1)[i],2,".png",sep=""), plot=ggps[[2]], width = 30, height = 30, units = "cm")
}
#ggps1 =.convertRatios(x,ratios_1,names(P0),main="Tumour Filtered to total",logx=logx, xlim = xlim)
#ggps2 =.convertRatios(x,ratios_2,names(P0),main="Unique to filtered",logx=logx, xlim = xlim)
#ggps3 =.convertRatios(x,ratios_3,names(P01),main="Benign filtered to total",logx=logx, xlim = xlim)
#ggps4 =.convertRatios(x,ratios_4,names(P0),main="Shared to total",logx=logx, xlim = xlim)

ggsave("ratio3.png", plot=ggps[[1]], width = 30, height = 30, units = "cm")


ggsave("ratio1.png", plot=ggps[[2]], width = 30, height = 30, units = "cm")

ncol1 = dim(P0)[2]
for(i in 1:ncol1){
  P0[,i] = P0[,i]/sum(P0[,i])
}

#density  = data.frame(P0)
#density = cbind(x,density[,!(names(density) %in% "combined")])
###DENSITY PLOTS



x = 2:1999
window=1
bins = NULL
#mergedT = readDataAll(names[1:4], '_P0.somatic_reads_read_length_count.json', x,"Tumour_filtered",window=window)
#s = apply(mergedT,1,sum)
#bins = .getBins(rev(s),x,5e2,rev=T)

P0 = readDataAll(names[1:4], '_P0_read_length_count.json', x,"Tumour_all",window=window, bins = bins)
shared = readDataAll(names[1:4], '_P0.shared.somatic_reads_read_length_count.json', x,"Tumour_shared",window=window, bins = bins) #'_P0_shared_read_length_count.json', x)
unique = readDataAll(names[1:4], '_P0.unique.somatic_reads_read_length_count.json', x, "Tumour_unique",window=window, bins = bins)
P01 = readDataAll(names[5:6], '_P0_read_length_count.json', x, "Benign_all",window=window, bins = bins)
merged = readDataAll(names[5:6], '_P0.somatic_reads_read_length_count.json', x, "Benign_filtered",window=window, bins = bins)
mergedT = readDataAll(names[1:4], '_P0.somatic_reads_read_length_count.json', x,"Tumour_filtered",window=window, bins = bins)
x = as.numeric(dimnames(P0)[[1]])

#vals = list()
#for(i in 1:dim(P0)[2]){
#  vals[[i]] = shared[[i]]/mergedT[[i]]
#}
#names(vals) = names(shared)
#df_vals = cbind(x,data.frame(vals))
#df2 = pivot_longer(df_vals, names_to = "ID", values_to = "ratio", cols=names(df_vals)[-1])
#ggp3<-ggplot(df2,aes(x, ratio,  fill=ID, color=ID,shape=ID))+ggtitle("Density plot(20bp sliding window)")+geom_line()

#a  = apply(shared,1,sum)
# b = apply(mergedT,1,sum)
# c  =a/b
# bins[which(c==max(c,na.rm=T)),]
# plot(x,c)
#secondpeak=33

smooth = 1
cumul = F
cumul1 =T
P0_d = data.frame(apply(P0, 2, makeDensity,smooth=smooth, cumul=cumul, cumul1=cumul1))
shared_d = data.frame(apply(shared, 2, makeDensity,smooth=smooth, cumul=cumul, cumul1=cumul1))
unique_d =data.frame( apply(unique, 2, makeDensity,smooth=smooth, cumul=cumul, cumul1=cumul1))
P01_d = data.frame(apply(P01, 2, makeDensity,smooth=smooth, cumul=cumul, cumul1=cumul1))
merged_d = data.frame(apply(merged,2,makeDensity,smooth=smooth, cumul=cumul, cumul1=cumul1))
merged_Td = data.frame(apply(mergedT,2,makeDensity,smooth=smooth, cumul=cumul, cumul1=cumul1))

#.addAvg(P0_d, "", probs  =probs)

l = list(Tumour_all =P0_d,Tumour_uniq= unique_d,Tumour_shared =shared_d,
         Benign_all=P01_d, Benign_filtered=merged_d,Tumour_filtered= merged_Td )
l = lapply(l, .addAvg,"",probs=probs)
.getDens<-function(x, l, nme,id_nme="ID"){
 res = list(x=x)
 for(i in 1:length(l)){
   res[[i+1]] = l[[i]][,which(names(l[[i]]) %in% nme)]
 }
 names(res) = c("x",names(l))
 density=data.frame(res)
 df2 = pivot_longer(density, names_to = id_nme, values_to = nme, cols=names(density)[-1])
 df2
}
inds = 1:length(l)  ## can change inds to include subset
#inds = 1:4
#inds = c(1,5)
#inds = 1:3
density = .getDens(x,l[inds],"X_q_50")
q_25 = .getDens(x,l[inds],"X_q_25", "ID")
q_50 = .getDens(x,l[inds],"X_q_50","ID1")
q_75 = .getDens(x,l[inds],"X_q_75","ID2")

df2 = cbind(q_25, q_50, q_75)[,c(1,2,3,6,9)]
df2$ID =as.factor(df2$ID)
na_ind = is.na(df2$X_q_50) | df2$X_q_50==0
df2 = df2[!na_ind,]

ggp2<-ggplot(df2,aes(x, X_q_50,  fill=ID, color=ID,shape=ID))+ggtitle("Density plot(20bp sliding window)")
ggp2<-ggp2+scale_y_continuous(trans='log10')
ggp2<-ggp2+geom_line()
#ggp2<-ggp2+xlim(0,200)
ggp2<-ggp2+theme_bw()+theme(text = element_text(size=textsize))+xlab("Fragment length(bp)")+ylab("Density")

#ggp2<-ggp2+geom_ribbon(aes(x,ymin=X_q_25, ymax=X_q_75, fill = ID),  alpha=0.1)
ggsave("density1.png", plot=ggp2, width = 30, height = 30, units = "cm")


