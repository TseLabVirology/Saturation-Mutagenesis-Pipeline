library("ggplot2")

data<- read.table("Sample.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)


sum=sum(data$Count)
data$Fraction <- data$Count/sum * 100

top<- data[1:1,]
top_sum=sum(top$Fraction)

other_sum=100-top_sum
otherdata<-data.frame("Other",1,0,other_sum)
names(otherdata)<-c("Sequence","Count","Rand","Fraction")

newdata <- rbind(top_sum,otherdata)

tiff("pie.tiff", units="in", width=5, height=5, res=300)

ggplot(newdata, aes(x="", y=Fraction, fill=Sequence))+ geom_bar(width = 2, stat = "identity", color="white") + coord_polar("y", start=0) + 
  theme_void() + scale_fill_manual(values=c("deeppink","black"))

dev.off()