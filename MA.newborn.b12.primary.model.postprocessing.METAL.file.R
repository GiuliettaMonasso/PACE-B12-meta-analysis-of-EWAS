infile<-"MA.b12.newborn.main.model.GSM.1.txt"
filepref<-"MA.b12.newborn.main.model.GSM"

#maken nieuwe df zonder HetDf = 0 
df<-read.table(infile, header = T, strings = F, dec=".", sep="\t")
summary(df)
 
df<-df[which(df$HetDf > 0),] 
summary(df)

# write file without row names, without quotes, and as csv
# outfile <- paste(filepref, ".csv", sep = "")
# write.csv(df, file = outfile, row.names = F, quote = F)

# read crossreactive probes: combination of Naeem and Chen
list<-read.csv("crossreactiveprobes.csv")
summary(list)
selection = subset(df, !(df$MarkerName %in% list$MarkerName))
summary(selection)
## Add Benjamini-Hochberg p-value correction to data

library(stats)
selection$P.value.FDR <- p.adjust(selection$P.value, method = "BH")
selection$P.value.FDR<-as.numeric(selection$P.value.FDR)
summary(selection)
head(selection)

# write.csv(selection, file = "MA.b12.newborn.main.model.GSM.1.zonderDf=0.cross.reactive.probes.csv", row.names = F, quote = F) 

# MA<-read.csv("MA.b12.newborn.main.model.GSM.1.zonderDf=0.cross.reactive.probes.csv", header = T, strings = F)

annotation<-"/home/860002/Methylation/Annotation/450k_annotationfile_v1_2.csv"

annotation<-read.csv(annotation, header = T, strings = F)
summary(annotation)
dim(annotation)
annotation<-annotation[,c(1:2,11:13,17:19,22:33)]
colnames(annotation)[2] <- "MarkerName"
summary(annotation)


MA<-merge(selection, annotation, by.x = "MarkerName", by.y = "MarkerName", all.x = T, all.y = F)
summary(MA)

MA<-MA[which(MA$CHR!="X"&MA$CHR!="Y"),]
summary(MA)




flagged<-"flaggedprobes.csv"
flagged<-read.csv(flagged, header = T, strings = F)
colnames(flagged)[1] <- "MarkerName"

MA<-merge(MA, flagged, by.x = "MarkerName", by.y = "MarkerName", all.x = T, all.y = F)
summary(MA)



outfile<- paste(filepref, ".csv", sep = "")
write.csv(MA, file = outfile, row.names = F, quote = F)

### make manhattan plot

Chr <- c(1:22,"X","Y")
GWplot<-function(data,P,Chr,title){
	print(summary(data))
        print(table(data$CHR))
        par(mar=c(5,5,4,2));
        phy.max<-tapply(data$MAPINFO, data$CHR,max,na.rm=T)
        cumlen=0
        for(i in Chr){

        data[data$CHR==i,"loc"]<-data[data$CHR==i,"MAPINFO"]+cumlen
        cumlen<-cumlen+phy.max[i]
        }
        phy.med<-tapply(data$loc, data$CHR,median,na.rm=T)[Chr]
        print(phy.med)
        data$mlgP<--log(data[,P], base=10)
        plot(data[,"loc"],data[,"mlgP"],type="n", xaxt="n", xlab="Chromosome",
            ylab="-Log p-value" ,main=title, cex.axis=1.5, cex.lab=2, xlim=c(0,max(data$loc,na.rm=T)))
        axis(side=1, at=phy.med[c(1:19)], labels=Chr[c(1:19)],
         tick=T,cex.axis=1.7,las=3)
        axis(side=1, at=phy.med[c(20:25)], labels=Chr[c(20:25)],
         tick=T,cex.axis=1.3,las=3,)
        for(i in Chr){
          if(which(Chr==i) %in% seq(2,30,2)) col="blue" else col="red"
          points(data[data$CHR==i,"loc"],data[data$CHR==i,"mlgP"],col=col,pch=20,cex=0.5)
          abline(7,0,col="grey", lty=2)
          abline(5.85,0,col="red", lty=3)
        }
}



observed <- sort(MA$P.value)
observed<- as.numeric(observed)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

outfile <- paste(filepref, ".QQplot.bmp", sep = "") 
bitmap(outfile, w=16,h=10)
par(cex = 1.3, cex.axis = 1.3, cex.lab = 1.15 )
plot(c(0,30), c(0,30), col="red", lwd=4, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,8), ylim=c(0,8), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black")
dev.off()


##  filter SE = 0 out to calculate lambda
## you could also have seen this in summary above! 

MA<-MA[which(MA$StdErr > 0),] 


MA$Z <- MA$Effect / MA$StdErr
MA$chi2 <- MA$Z * MA$Z
lambda<-median(MA$chi2)/0.455

outfile<- paste(filepref,".lambda",".txt")
write.table(lambda,file=outfile, row.names=F, quote=F)

## checked for probes with information from only 1 study (zero degrees of freedom)

# MA2<-MA[which(MA$HetDF == 0),]
# dim (MA2)


outfile<- paste(filepref, "logplot.bmp", sep = "")
bitmap(outfile, w=16,h=10)
GWplot(MA,"P.value",Chr,title="")
dev.off()










