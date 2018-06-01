rm(list=ls())
data<-read.csv('~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/Exons.GeneticDistance.LogBins.csv')
data$B<-data$pi/0.0083
data<-data[order(abs(data$mid)),]

loessMod10 <- loess(B ~ abs(mid), data=data, span=0.2, weights = sites) # 10% smoothing span
smoothed10 <- predict(loessMod10) 

plot(data$B, x=abs(data$mid), type="p", log='x',xlim=c(1,1500))
lines(smoothed10, x=abs(data$mid), col="red")
data$Bsmooth <- smoothed10

write.csv(data, file = '~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/Exons.GeneticDistance.LogBins.Loess.csv')

#########
#########
#########

rm(list=ls())

data<-read.csv('~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/CNEs.GeneticDistance.LogBins.csv')
data$B<-data$pi/0.0083
data<-data[order(abs(data$mid)),]

loessMod10 <- loess(B ~ abs(mid), data=data, span=0.2, weights = sites) # 10% smoothing span
smoothed10 <- predict(loessMod10) 

plot(data$B, x=abs(data$mid), type="p", xlim=c(1,250), log = 'x')
lines(smoothed10, x=abs(data$mid), col="red")
data$Bsmooth <- smoothed10

write.csv(data, file = '~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/CNEs.GeneticDistance.LogBins.Loess.csv')

rm(list=ls())

exon <- read.csv( '~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/Exons.GeneticDistance.LogBins.Loess.csv')
exon$label <- 'Protein-Coding Exon'
exon <- exon[abs(exon$mid) < 2500, ]
cne <- read.csv( '~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/CNEs.GeneticDistance.LogBins.Loess.csv')
cne$label <- 'Conserved Non-Coding Element'
cne <- cne[abs(cne$mid) < 250, ]

da <- rbind(exon, cne)
da <- da[abs(da$mid) > 1, ]
library(ggplot2)

ggplot(data = da, aes(x= abs(mid), y = B))+
  geom_point(lwd = 0.8)+
  geom_line(aes(x = abs(mid), y= Bsmooth), col = 'red', lwd = 0.8)+
  facet_grid(~label, scales = 'free_x')+
  scale_x_log10()+
  xlab('Distance from Element (4Ner)')+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=0,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.x = element_text(size = 15)
  )

