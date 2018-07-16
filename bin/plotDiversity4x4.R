rm(list=ls())
library(ggplot2)
cbPalette <- c("#56B4E9",  "#D55E00", "#E69F00")


cast_noGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/exons/autosomes.ncpg.noGC.castaneus.csv')
cast_noGC$label<-'No Gene Conversion'
cast_pGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/exons/autosomes.ncpg.PaigenGC.castaneus.csv')
cast_pGC$label<-'Low Gene Conversion'
cast_hGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/exons/autosomes.ncpg.HighGC.castaneus.csv')
cast_hGC$label<-'High Gene Conversion'
z<-rbind(cast_noGC, cast_pGC, cast_hGC)
z$map <- 'LD-based'

Cox_noGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/exons/autosomes.ncpg.noGC.Cox.csv')
Cox_noGC$label<-'No Gene Conversion'
Cox_pGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/exons/autosomes.ncpg.PaigenGC.Cox.csv')
Cox_pGC$label<-'Low Gene Conversion'
Cox_hGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/exons/autosomes.ncpg.HighGC.Cox.csv')
Cox_hGC$label<-'High Gene Conversion'
w<-rbind(Cox_noGC, Cox_pGC, Cox_hGC)
w$map <- 'Pedigree-based'

z<-rbind(z,w)
z$element<-'Protein-Coding Exons'
z<-z[z$mid < 2500,]
z<-z[z$mid > 1,]

cast_noGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/cneAnalysis/autosomes.ncpg.NoGC.castaneus.csv')
cast_noGC$label<-'No Gene Conversion'
cast_pGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/cneAnalysis/autosomes.ncpg.PaigenGC.castaneus.csv')
cast_pGC$label<-'Low Gene Conversion'
cast_hGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/cneAnalysis/autosomes.ncpg.HighGC.castaneus.csv')
cast_hGC$label<-'High Gene Conversion'
f<-rbind(cast_noGC, cast_pGC, cast_hGC)
f$map <- 'LD-based'

Cox_noGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/cneAnalysis/autosomes.ncpg.NoGC.Cox.csv')
Cox_noGC$label<-'No Gene Conversion'
Cox_pGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/cneAnalysis/autosomes.ncpg.PaigenGC.Cox.csv')
Cox_pGC$label<-'Low Gene Conversion'
Cox_hGC<-read.csv('PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/both/cneAnalysis/autosomes.ncpg.HighGC.Cox.csv')
Cox_hGC$label<-'High Gene Conversion'
g<-rbind(Cox_noGC, Cox_pGC, Cox_hGC)
g$map <- 'Pedigree-based'

q<-rbind(f,g)
q$element<-'Conserved Non-Coding Elements'
q<-q[q$mid < 250,]
q<-q[q$mid > 1,]

v <- rbind(q,z)
ggplot(data = v, aes(x = mid, y = pi, col = label))+
  geom_line(lwd = 1, alpha = 0.8)+
 #   scale_x_log10(limits = c(1,300))+
  scale_y_continuous(expression(pi))+
  scale_x_continuous(expression('Distance From Element (4'*N[e]*'r)'))+
  facet_grid(map~element, scales= 'free_x')+
  scale_color_manual('',values = cbPalette)+
  theme_bw()+
  theme(
    #    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0, face = 'italic'),
    axis.title.y = element_text(size=20,angle=0,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.x = element_text(size = 13, face = 'italic'),
    strip.text.y = element_text(size = 15)
  )

ggplot(data = v, aes(x = mid, y = rat_div_jc, col = label))+
  geom_line(lwd = 1)+
  #   scale_x_log10(limits = c(1,300))+
  scale_y_continuous(expression(pi), limits = c(0.1,.2))+
  # scale_x_continuous(limits = c(1,300))+
  facet_grid(map~element, scales= 'free_x')+
  scale_color_manual('',values = cbPalette)+
  theme_bw()+
  theme(
    #    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0, face = 'italic'),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.x = element_text(size = 15, face = 'italic'),
    strip.text.y = element_text(size = 15)
  )

ggplot(data = v, aes(x = mid, y = pi/rat_div_jc, col = label))+
  geom_line(lwd = 1)+
  scale_y_continuous(expression(pi))+
  scale_x_continuous()+
  facet_grid(map~element, scales= 'free_x')+
  scale_color_manual('',values = cbPalette)+
  theme_bw()+
  theme(
    #    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0, face = 'italic'),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.x = element_text(size = 15, face = 'italic'),
    strip.text.y = element_text(size = 15)
  )

ggplot(data = v, aes(x = mid, y = sites/1e6, col = label))+
  geom_line(lwd = 1)+
  scale_y_log10('Sites (Mbp)')+
  scale_x_continuous()+
  facet_grid(map~element, scales= 'free_x')+
  scale_color_manual('',values = cbPalette)+
  theme_bw()+
  theme(
    #    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0, face = 'italic'),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.x = element_text(size = 15, face = 'italic'),
    strip.text.y = element_text(size = 15)
  )
