rm(list=ls())


library(ggplot2)

exon<-read.csv('~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/Exons.Genetic.Logbins.csv')
exon$B<-exon$pi/0.0083
exon<-exon[order(abs(exon$mid)),]

loessMod10.exon <- loess(B ~ abs(mid), data=exon, span=0.2, weights = sites) # 10% smoothing span
smoothed10.exon <- predict(loessMod10.exon, se=T) 

exonD <- data.frame(cbind(smoothed10.exon$fit, smoothed10.exon$s, abs(exon$mid), abs(exon$B)))
exonD$lab<- 'Protein-Coding Exons'
exon$Bsmooth <- smoothed10.exon$fit

write.csv(exon, file = '~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/Exons.GeneticDistance.LogBins.Loess.csv')

#########
#########
#########


cne<-read.csv('~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/CNEs.GeneticDistance.LogBins.csv')
cne$B<-cne$pi/0.0083
cne<-cne[order(abs(cne$mid)),]

loessMod10.cne <- loess(B ~ abs(mid), data= cne, span=0.2, weights = sites) # 10% smoothing span
smoothed10.cne <- predict(loessMod10.cne, se= T) 

plot(cne$B, x=abs(cne$mid), type="p", xlim=c(1,250), log = 'x')
lines(smoothed10$fit, x=abs(cne$mid), col="red")

cneD <- data.frame(cbind(as.numeric(smoothed10.cne$fit), as.numeric(smoothed10.cne$s), as.numeric(abs(cne$mid)), as.numeric(abs(cne$B))))
cneD$lab <- 'Conserved Non-Coding Elements'

cne$Bsmooth <- smoothed10.cne$fit
write.csv(cne, file = '~/PhD/Coding/Estimate_selection_from_troughs/BGS/analysis/CNEs.GeneticDistance.LogBins.Loess.csv')
str(cneD)

cairo_pdf('~/PhD/Coding/Estimate_selection_from_troughs/BGS/BGSplotLoess.pdf', height = 8, width = 10)

ggplot(data = cneD, aes(x= X3, y = X4))+
  geom_point(lwd = 0.8)+
  geom_line(aes(x= X3, y = X1), col = 'red')+
  geom_ribbon(aes(x= X3, ymax = X1+(2*X2), ymin = X1-(2*X2)), fill = 'red', alpha = 0.35)+
  geom_point(data = exonD, aes(x= X3, y = X4),lwd = 0.8)+
  geom_line(data = exonD, aes(x= X3, y = X1), col = 'red')+
  geom_ribbon(data = exonD, aes(x= X3, ymax = X1+(2*X2), ymin = X1-(2*X2)), fill = 'red', alpha = 0.35)+
  facet_grid( ~ lab, scales = 'free_x')+
  xlab(expression('Distance from Element (4'*N[e]*'r)'))+
  ylab('B')+
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
dev.off()
  

  