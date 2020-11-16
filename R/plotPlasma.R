.libPaths("C:/Users/LCOIN/R-4.0.2/library")

library('jsonlite');
library('ggplot2')
library(tidyr)
library(reshape2)
RHOME="../../R"
setwd("../data/new_samples")
source(paste(RHOME, "plotPlasmaFuncs.R", sep="/"));

names = c("1524", "1249", "1494", "1084","065", "098");




#x = 2:1999
x = 2:1999
P0 = readDataAll(names[1:4], '_P0_read_length_count.json', x)
shared = readDataAll(names[1:4], '_P0.shared.somatic_reads_read_length_count.json', x) #'_P0_shared_read_length_count.json', x)
unique = readDataAll(names[1:4], '_P0.unique.somatic_reads_read_length_count.json', x)

.cumsum1<-function(x) sum(x) -cumsum(x)

P0_c = data.frame(cbind(apply(P0,2,cumsum), apply(P0,2,.cumsum1)))
shared_c = data.frame(cbind(apply(shared,2,cumsum), apply(shared,2,.cumsum1)))
unique_c = data.frame(cbind(apply(unique,2,cumsum), apply(unique,2,.cumsum1)))
ncol = dim(P0_c)[2]
ratios = matrix(nrow = dim(P0_c)[1],ncol = ncol)
for(i in 1:ncol){

  shared_c[,i] = shared_c[,i]/P0_c[,i]
  unique_c[,i] = unique_c[,i]/P0_c[,i]
  ratios[,i] = shared_c[,i]/unique_c[,i]
}

ratios = data.frame(ratios)
names(ratios) = c(paste(names(unique),"less than or equal to",sep="."),paste(names(unique),"greater than",sep="."))
ratios = cbind(x,ratios)
df = pivot_longer(ratios, names_to = "ID", values_to = "proportion", cols=names(ratios)[-1])  %>% 
      separate(ID, c('ID', 'type'), sep='\\.', remove = T)
ratios1 = ratios[,c(1,grep("combined",names(ratios)))]
df1 = pivot_longer(ratios1, names_to="ID",values_to="proportion",cols = names(ratios1)[-1])
ggp<-ggplot(df1, aes(x,proportion, linetype=ID))+geom_line()+ggtitle("Ratio of shared to unique")+ylab("Ratio")
ggp<-ggp+scale_x_continuous(trans='log10')
textsize=20
ggp<-ggp+theme_bw()+theme(text = element_text(size=textsize))
ggsave("ratio.png", plot=ggp, width = 30, height = 30, units = "cm")

ggp1 = ggplot(df, aes(x,proportion, fill=ID,color=ID, linetype=type))+geom_line()+scale_x_continuous(trans='log10')+ggtitle("Ratio of shared to unique")
ggp1<-ggp1+theme_bw()+theme(text = element_text(size=textsize))
ggsave("ratio1.png", plot=ggp1, width = 30, height = 30, units = "cm")

ncol1 = dim(P0)[2]
for(i in 1:ncol1){
  P0[,i] = P0[,i]/sum(P0[,i])
}

#density  = data.frame(P0)
#density = cbind(x,density[,!(names(density) %in% "combined")])
###DENSITY PLOTS

x = 2:1999
P0 = readDataAll(names[1:4], '_P0_read_length_count.json', x)
shared = readDataAll(names[1:4], '_P0.shared.somatic_reads_read_length_count.json', x) #'_P0_shared_read_length_count.json', x)
unique = readDataAll(names[1:4], '_P0.unique.somatic_reads_read_length_count.json', x)
P01 = readDataAll(names[5:6], '_P0_read_length_count.json', x)
merged = readDataAll(names[5:6], '_P0.somatic_reads_read_length_count.json', x)
mergedT = readDataAll(names[1:4], '_P0.somatic_reads_read_length_count.json', x)


P0_d = data.frame(apply(P0, 2, makeDensity,smooth=10))
shared_d = data.frame(apply(shared, 2, makeDensity,smooth=10))
unique_d =data.frame( apply(unique, 2, makeDensity,smooth=10))
P01_d = data.frame(apply(P01, 2, makeDensity,smooth=10))
merged_d = data.frame(apply(merged,2,makeDensity,smooth=10))
merged_Td = data.frame(apply(merged,2,makeDensity,smooth=10))



density = data.frame(x=x,"Tumour_all" = P0_d$combined,
          "Tumour_uniq" = unique_d$combined, 
          "Tumour_shared" = shared_d$combined, 
          "Benign_all" = P01_d$combined,
          "Benign_filtered" = merged_d$combined, 
          "Tumour_filtered" = merged_Td$combined)


df2 = pivot_longer(density, names_to = "ID", values_to = "proportion", cols=names(density)[-1])

ggp2<-ggplot(df2, aes(x, proportion, fill=ID, color=ID))+geom_line()+scale_y_continuous(trans='log10')+ggtitle("Density plot(20bp sliding window)")
ggp2<-ggp2+theme_bw()+theme(text = element_text(size=textsize))+xlab("Fragment length(bp)")+ylab("Density")
ggsave("density.png", plot=ggp2, width = 30, height = 30, units = "cm")


###FROM HERE IS OLDER SCRIPT



#li = list(unique_d = unique_d[,4,drop=F],shared_d = shared_d[,4,drop=F])
resdir = "res"
dir.create(resdir)
pdf(paste(resdir,"out.pdf",sep="/"))
plotAll1(x, li, n=1, log = "", inds = 1:500, func=cumsum)
plotAll1(x, li, n=1, log = "y", inds = 1:1998, func=cumsum1)
plotAll1(x, li, n=20, log = "", inds = 1:1998, func=ident)
plotAll1(x, li, n=20, log = "x", inds = 1:1998, func=ident)

plotAll1(x, li, n=20, log = "y", inds = 1:1998, func=ident)
plotAll1(x, li, n=20, log = "xy", inds = 1:1998, func=ident)
dev.off()




thresh = 1:2000
secondThresh = thresh - 10;
l = plotAll(outf= NULL, thresh = thresh, lt =F , secondThresh = secondThresh, ylim = c(1e-5,1))

replot(l, xlim = c(2,1998), norm=F, smooth = 50)


replot(l, xlim = c(1,2000), norm=F, smooth = 10)
l = plotAll(outf= paste(resdir,"plots.pdf",sep="/"), thresh = thresh, lt=T, secondThresh = 0, ylim = c(1e-5,1))
thresh = 1:200
secondThresh = thresh - 10;
l = plotAll(outf=paste(resdir,"plots1.pdf",sep="/"), lt=T, thresh=thresh, secondThresh = secondThresh)

l = plotAll(outf=NULL, lt=T, secondThresh = 20, thresh=1:200, norm=NULL, ylim = c(0,5))
a1 = plotlines(average(l), outf = paste(resdir,"average_short_20.pdf",sep="/"))
plotDiff(a1, l[[1]]$thresh, outf=paste(resdir,"average_diff_short_20.pdf",sep="/"))


#a1 = plotlines(l)



lGT = plotAll(outf=paste(resdir,'greaterThan.pdf',sep="/"), lt=F, secondThresh = 2000, thresh=0:2000, norm=NULL, ylim = c(0.00001,0.05))
lLT = plotAll(outf=paste(resdir, 'lessThan.pdf',sep="/"), lt=T, secondThresh = 0, thresh=0:300, norm=NULL, ylim = c(0.00001,0.05))
a1 = plotlines(average(l), outf = paste(resdir, "average_long.pdf",sep="/"))
plotDiff(a1, l[[1]]$thresh, outf=paste(resdir,"average_diff_long.pdf",sep="/"))


a1 = plotlines(l)

norm = plotAll(T, outf=paste(resdir,"normalised.pdf",sep="/"), thresh = 1:200)
#


l100 = plotAll(outf=NULL, lt=T, secondThresh = 20, thresh=1:200, norm=NULL, ylim = c(0,5))
l50 = plotAll(outf=NULL, lt=T, secondThresh = 20, thresh=1:200, norm=NULL, ylim = c(0,5), targetinds = c(50,200))


getC(l$c1524s$w,l$c1524$rates, mutrates = c(0.1, 1e-5))
getC(l$c1524s$w,l$c1524s$rates, mutrates = c(0.1, 1e-5))

getC(l$c1249s$w,l$c1249$rates, mutrates = c(0.1, 1e-5))
 getC(l$c1494s$w,l$c1494$rates, mutrates = c(0.1, 1e-5))


w = apply(l$c1524$c[c(100,200),],2,normEnd)[1,]

#w[2] = apply(l$c616$c[c(100,200),],2, normEnd)[1,1]  ##perc fragments in normal plasma lt 100bp
#w[1] = apply(l$c1494$c[c(100,200),],2, normEnd)[1,2]  ##perc fragments in tumour plasma lt 100bp


calcC(w,ratio(apply(l$c1524$c[c(100,200),],1,ratio)))


100*calcC(w,ratio(apply(l$c1524$c[c(100,200),],1,ratio)))
100*calcC(w,ratio(apply(l$c1494$c[c(100,200),],1,ratio)))
100*calcC(w,ratio(apply(l$c1249$c[c(100,200),],1,ratio)))
100*calcC(w,ratio(apply(l$c606$c[c(100,200),],1,ratio)))
100*calcC(w,ratio(apply(l$c616$c[c(100,200),],1,ratio)))


#a = calcDiff(c1249)
#a = calcDiff(c1524)
#a = calcDiff(c1494)
# rat1524 = apply(c1524,1,ratio)
# rat1494 = apply(c1494,1,ratio)
# rat1249 = apply(c1249,1,ratio)
#((( rat[2]) * 0.07)* (1/rat[1]))

c606 = readData(ext_all='616_P0_read_length_count.json', ext_unique = '609_P0_read_length_count.json', thresh = thresh )






#####
.makeCumlative<-function(mat){
 mat1 = data.frame(apply(mat,2,cumsum))
}



