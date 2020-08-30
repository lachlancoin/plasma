library('jsonlite');
RHOME="../R"
source(paste(RHOME, "plotPlasmaFuncs.R", sep="/"));


names = c("1524", "1249", "1494", "609", "616");




#x = 2:1999
x = 2:1999
P0 = readDataAll(names[1:3], '_P0_read_length_count.json', x)
shared = readDataAll(names[1:3], '_P0_shared_read_length_count.json', x)
unique = readDataAll(names[1:3], '_P0_unique_read_length_count.json', x)
mergedT = readDataAll(names[1:3], '_P0_merged_read_length_count.json', x)
P01 = readDataAll(names[4:5], '_P0_read_length_count.json', x)
merged = readDataAll(names[4:5], '_P0_merged_read_length_count.json', x)
P0_d = apply(P0, 2, makeDensity)
shared_d = apply(shared, 2, makeDensity)
unique_d = apply(unique, 2, makeDensity)
P01_d = apply(P01, 2, makeDensity)
merged_d = apply(merged,2,makeDensity)
merged_Td = apply(mergedT,2,makeDensity)

li = list("Tumour_all" = P0_d[,4,drop=F], "Tumour_uniq" = unique_d[,4,drop=F], "Tumour_shared" = shared_d[,4,drop=F], "Benign_all" = P01_d[,3,drop=F], "Benign_filtered" = merged_d[,3,drop=F], "Tumour_filtered" = merged_Td[,4,drop=F])

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



