#Libraries
library(reshape2)
library(ggplot2)
library(vegan)

#Graphs
phylum<-read.csv(file="Data/summer19_noEuk_OTU.csv", header=T)
mapping<-read.csv(file="Data/mapping_summer19.csv")

head(phylum)
head(mapping)

p.m<-melt(phylum)

dates<-as.vector(mapping$Date)

#Phylum Barplot
ggplot(p.m,aes(x=variable, y=value, fill=Phylum)) +geom_bar(stat = "identity") +
  xlab("Sample") +
  ylab("Relative Abundances (%)") +
  theme_bw(base_size = 20) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.text = element_text(angle = 90)) + scale_fill_manual(values = c("#FCD16B","#C7CEF6","#D3DDDC","#C6B19D","#B62A3D","#B5966D","#CECD7B","#F24D29","#A35E60","#CC8B3C","#957A6D","#AC6E49","#9A872D","#1C366B","#D1362F","#456355","#DBB165","#FDDDA4", "#C7CEF6","#27223C"))

#Unifrac PcoA (Unweighted)
unifrac<-read.csv(file="Data/undistance-matrix.csv", header=T, row.names = 1)
unifrac.ord<-pcoa(unifrac)
uniord.s<-unifrac.ord$vectors
write.csv(uniord.s, file="Data/scores2.csv")
uniord.s<-read.csv(file="Data/scores2.csv", check.names = T)

ggplot(data.frame(uniord.s),aes(x=Axis.1, y=Axis.2, shape=Air.Type, color=Type)) + geom_point(size=7) +
  xlab("PcoA 1 (26.20%)") +
  ylab("PcoA 2 (14.90% )") +
  theme_bw(base_size = 20) +
  scale_color_manual(values=c("#F24D29", "#1DACE8", "#1C366B")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"   ))+theme(axis.text = element_text(angle = 90)) 

ps.disper <- betadisper(as.dist(unifrac), mapping$Type)
permutest(ps.disper, pair=TRUE)

# Function PcoA

fd<-read.csv(file="Data/pcoa_func.csv", header=T) #PCOA FROM QIIME2 PICRUST

ggplot(data.frame(fd),aes(x=Axis.1, y=Axis.2, shape=Air.Type, color=Type)) + geom_point(size=7) +
  xlab("PcoA 1 (65.44%)") +
  ylab("PcoA 2 (19.66% )") +
  theme_bw(base_size = 20) +
  scale_color_manual(values=c("#76A08A", "#541F12", "#A35E60")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  theme(axis.text = element_text(angle = 90)) 


