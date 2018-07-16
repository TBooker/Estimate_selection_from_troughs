rm(list=ls())

library(ggplot2)
library(wesanderson)
options(scipen=10000)
TommyTheme <-   
  theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=17,vjust=0.5, angle = 0),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    strip.text.y = element_text(size = 20, face = 'italic'),
    strip.text.x = element_text(size = 20, face = 'italic')
  )


pal <- wes_palette(n = 5, name = "Royal1", type = "continuous")

ex<-read.csv('/Users/s0784966//PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/bestExonModel.csv')
ex$element<-'Protein-Coding Exons'
cn<-read.csv('/Users/s0784966//PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/bestCNEModel.csv')
cn$element<-'Conserved Non-Coding Elements'

z<- rbind(ex,cn)

#z$element <- factor(z$element, levels = c('Protein-Coding Exons', 'Conserved Non-Coding Elements'), labels = c('Protein-Coding Exons', 'Conserved Non-Coding Elements'))
#
cairo_pdf('/Users/s0784966//PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/ModelFit.pdf', width = 12, height = 6)

ggplot(data = z, aes(x = distance, y = pi, col = 'Observed'))+
  geom_line(lwd = 1.1)+
  geom_line(aes(x = distance, y = BGS, col = 'B'), lwd = 1.05, lty = 'dashed')+
  geom_line(aes(x = distance, y = fitted, col = 'Fitted'), lwd = 1.2, lty = 'dashed')+
  xlab(expression('Distance to Element (4'*N[e]*'r)'))+
  facet_grid(~element,scales = 'free_x')+
  scale_color_manual('',values = pal)+
  ylab(expression(pi/pi[0]))+
  theme_bw()+
  TommyTheme

dev.off()
