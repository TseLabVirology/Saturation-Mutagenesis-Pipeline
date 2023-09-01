library("ggplot2")

data<- read.table("merge.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

sel <- subset(data, Rand > 100)
unsel <- subset(data, Rand < 100)
#pull out only unselected values, ordered by sel counts

sel_sum=sum(sel$Count)
unsel_sum=sum(unsel$Count)

sel$Fraction <- sel$Count/sel_sum * 100
unsel$Fraction <- unsel$Count/unsel_sum *100

sel$Enriched <- sel$Fraction/unsel$Fraction
unsel$Enriched <- unsel$Fraction

newdata <-rbind(sel,unsel)

sel_sort <- sel[order(-sel$Enriched),]
#order selected values by enrichment
top_sel_sort <-sel_sort[1:5,]
unsel_sort <- unsel[order(-unsel$Count),]
#order unselected values by their counts
top_sel<- data[1:3,]
top_unsel <- unsel[1:3,]
top_sort <- unsel_sort[1:1,]
scale=c(10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2))
label=c(-5,-4,-3,-2,-1,0,1,2)
windowsFonts(A=windowsFont("Arial"))
c<-qplot(x=Rand, y=Fraction, data=newdata, color=I("deeppink"), size=Enriched, alpha=0.5, geom="point", xlab="Fusion Loop Library", ylab="log(Fraction of Reads)")

tiff("Sample.tiff", units="in", width=10, height=8, res=300)
c + scale_radius(range=c(1,20)) + scale_x_continuous(breaks=c(50,150), labels=c("Unselected","Selected")) + 
  scale_y_continuous(trans = "log", limits=c(0.00001,100), breaks=scale, labels=label) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_text(family="A", color="black", size=20),
        axis.text.y=element_text(family="A", color="black", size=20), axis.title.x=element_text(family="A", color="black",size=24), 
        axis.title.y=element_text(family="A", color="black",size=24), legend.position="none", 
        panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),panel.background=element_rect(fill="grey98")) + 
  #annotate("text", x=as.numeric(top_unsel$Rand), y=as.numeric(top_unsel$Fraction), family="A", label=top_unsel$Sequence ) + 
  #annotate("text", x=as.numeric(top_sel_sort$Rand), y=as.numeric(top_sel_sort$Fraction), family="A", label=top_sel_sort$Sequence) + 
  geom_vline(xintercept=100)
dev.off()
