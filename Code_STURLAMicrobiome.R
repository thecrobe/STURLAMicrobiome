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


#Sturla Class Distributions
sturla<-read.csv(file="Data/STURLA_Mean_Class.csv", header=T)
sturla$Class<-factor(sturla$Class, levels = sturla$Class)

ggplot(sturla, aes(x=sturla$Class, y=sturla$Mean, fill=sturla$Class)) +geom_bar(stat="identity")  + theme_bw(base_size = 20) + scale_color_manual(values=c("#A0D694", "#9359DB", "#CC6553", "#475E7F", "#EDA4DD","#D8B461","")) +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +theme(axis.text = element_text(angle = 90)) 

#STURLA Correlations
data<-read.csv("Data/STURLA_CORR.csv", header=T,row.names = 1)
cor(data,data,method="spearman") # correlation matrix

library(RVAideMemoire)
perm.cor.test(data$Function_Shannon,data$tgpl)
perm.cor.test(data$Function_Shannon,data$tgplm)
perm.cor.test(data$Function_Shannon,data$tgp)
perm.cor.test(data$Function_Shannon,data$tgwp)
perm.cor.test(data$Function_Shannon,data$tgbp)
perm.cor.test(data$Function_Shannon,data$tgbpl)
perm.cor.test(data$Function_Shannon,data$tgplmh)
perm.cor.test(data$Function_Shannon,data$tgwpl)
perm.cor.test(data$Function_Shannon,data$twpm)
perm.cor.test(data$Function_Shannon,data$other)

perm.cor.test(data$Phylogenetic,data$tgpl)
perm.cor.test(data$Phylogenetic,data$tgplm)
perm.cor.test(data$Phylogenetic,data$tgp)
perm.cor.test(data$Phylogenetic,data$tgwp)
perm.cor.test(data$Phylogenetic,data$tgbp)
perm.cor.test(data$Phylogenetic,data$tgbpl)
perm.cor.test(data$Phylogenetic,data$tgplmh)
perm.cor.test(data$Phylogenetic,data$tgwpl)
perm.cor.test(data$Phylogenetic,data$twpm)
perm.cor.test(data$Phylogenetic,data$other)

perm.cor.test(data$Observed,data$tgpl)
perm.cor.test(data$Observed,data$tgplm)
perm.cor.test(data$Observed,data$tgp)
perm.cor.test(data$Observed,data$tgwp)
perm.cor.test(data$Observed,data$tgbp)
perm.cor.test(data$Observed,data$tgbpl)
perm.cor.test(data$Observed,data$tgplmh)
perm.cor.test(data$Observed,data$tgwpl)
perm.cor.test(data$Observed,data$twpm)
perm.cor.test(data$Observed,data$other)

perm.cor.test(data$Pielou.s.Evenness,data$tgpl)
perm.cor.test(data$Pielou.s.Evenness,data$tgplm)
perm.cor.test(data$Pielou.s.Evenness,data$tgp)
perm.cor.test(data$Pielou.s.Evenness,data$tgwp)
perm.cor.test(data$Pielou.s.Evenness,data$tgbp)
perm.cor.test(data$Pielou.s.Evenness,data$tgbpl)
perm.cor.test(data$Pielou.s.Evenness,data$tgplmh)
perm.cor.test(data$Pielou.s.Evenness,data$tgwpl)
perm.cor.test(data$Pielou.s.Evenness,data$twpm)
perm.cor.test(data$Pielou.s.Evenness,data$other)

cor.matrix<-read.csv(file = "Data/sturladivcorrelationrho.csv")
cor.melt<-melt(cor.matrix)

ggplot(cor.melt, aes(x=cor.melt$variable, cor.melt$Correlation, fill=cor.melt$value)) +geom_tile() + scale_fill_continuous(type = "viridis") +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 


#Mantel Test
sturla<-read.csv(file="./Data/sturla_full.csv", header=T, row.names = 1)
sturla.bray<-vegdist(t(sturla))

unifrac.nosub<-read.csv(file = "Data/undistance-matrix-nosubway.csv",header = T, row.names = 1)

dim(sturla.bray)
set.seed(3)
mantel(sturla.bray,unifrac.nosub)
