library("ggplot2")
library("dplyr")

data<- read.table("mergeAA.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

data$Fuzzy <- as.integer(as.numeric(data$Fuzzy))
data <- data %>% mutate(AA_bin = cut(Fuzzy, breaks=c(0,20,35,55,70,85,110), right=FALSE))
#depending on the number of amino acids in the sequence, the breaks will need to be changed. Fuzzy value corresponds to percentage similarity, adjust bins to correspond to 1 amino acid intervals.

unsel <- subset(data, Rand < 100) #pull out only unselected values, ordered by sel counts
P1 <- subset(data, Rand > 100)

P1_sum=sum(P1$Count)
unsel_sum=sum(unsel$Count)

P1$Fraction <- P1$Count/P1_sum * 100
unsel$Fraction <- unsel$Count/unsel_sum * 100

P1$Enriched <- P1$Fraction/unsel$Fraction
unsel$Enriched <- unsel$Fraction

P1_sort <- P1[order(-P1$Enriched),]
unsel_sort <- unsel[order(-unsel$Enriched),]

top_P1<- P1_sort[1:5,]
top_unsel <- unsel_sort[1:3,]

windowsFonts(A=windowsFont("Arial"))

tiff("bubble-AA.jpeg", units="in", width=5, height=8, res=300)

qplot(x=Rand, y=Count, data=P1, color=AA_bin, size=Enriched, alpha=0.5, geom="point", xlab= "X Axis Label", ylab="log(Reads)") +
  scale_radius(range=c(1,20)) + scale_x_continuous(breaks=c(150), labels=c("labels")) + 
  scale_y_continuous(trans = "log", limits=c(1,10000), breaks=c(10^0,10^1,10^2,10^3,10^4), labels=c(0,1,2,3,4)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_text(family="A", color="black", size=14), 
        axis.text.y=element_text(family="A", color="black", size=14), axis.title.x=element_text(family="A", color="black", size=20), 
        axis.title.y=element_text(family="A", color="black", size=20), legend.position="right", 
        panel.border=element_rect(colour="black", fill=NA, size=1),panel.background=element_rect(fill="grey98"))  +
  scale_color_manual(values = c("darkgray","lightpink","seagreen","gold","skyblue3","black")) +
  annotate("text", x=as.numeric(top_P1$Rand), y=as.numeric(top_P1$Count), family="A", label=top_P1$Sequence)

dev.off()

