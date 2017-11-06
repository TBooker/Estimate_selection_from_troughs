rm(list=ls())
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

all<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/All.Rho.csv')
all$source<-'All Data'
allBGS<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/BGS.All.Rho.csv')
allBGS$source<-'All Data - BGS only'

plotter <- rbind( all, allBGS )

ggplot(data = plotter, aes( x = dist/1000, y= pi, col = source))+
  geom_line()+
  ylab(expression(pi))+
  xlab('Distance to Exon (Kbp)')+
  scale_color_manual('', values = cbPalette)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17, angle = 0, vjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13)
  )


#######################################################################
cbPalette <- c("#56B4E9", "#009E73","#0072B2", "#D55E00", "#CC79A7")

r0.009<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/Rho0.009.csv')
r0.009$source<-'Rho = 0.009'
r0.009$rho <- r0.009$dist*0.009
r0.0045<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/Rho0.0045.csv')
r0.0045$source<-'Rho = 0.0045'
r0.0045$rho <- r0.0045$dist*0.0045
r0.0009<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/Rho0.0009.csv')
r0.0009$source<-'Rho = 0.0009'
r0.0009$rho <- r0.0009$dist*0.0009

plotter <- rbind( r0.009, r0.0045, r0.0009)


ggplot(data = plotter, aes( x = dist/1000, y= pi, col = source))+
  geom_line()+
  ylab(expression(pi))+
  xlab('Distance to Exon (Kbp)')+
  scale_color_manual('', values = cbPalette)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17, angle = 0, vjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

ggplot(data = plotter, aes( x = rho, y= pi, col = source))+
  geom_line()+
  ylab(expression(pi))+
  xlab('Distance to Exon (Rho)')+
  scale_color_manual('', values = cbPalette)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17, angle = 0, vjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

###################################################################################

br0.009<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/BGS.Rho0.009.csv')
br0.009$source<-'Rho = 0.009 - BGS Only'
br0.009$rho <- br0.009$dist*0.009
br0.0045<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/BGS.Rho0.0045.csv')
br0.0045$source<-'Rho = 0.0045 - BGS Only'
br0.0045$rho <- br0.0045$dist*0.0045
br0.0009<-read.csv('project/15.Estimate_Selection_From_Trough/simulations/Nes95_dDFE.Exponential/BGS.Rho0.0009.csv')
br0.0009$source<-'Rho = 0.0009 - BGS Only' 
br0.0009$rho <- br0.0009$dist*0.0009

plotter <- rbind( br0.009, br0.0045, br0.0009)


ggplot(data = plotter, aes( x = dist/1000, y= pi, col = source))+
  geom_line()+
  ylab(expression(pi))+
  xlab('Distance to Exon (Kbp)')+
  scale_color_manual('', values = cbPalette)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17, angle = 0, vjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

ggplot(data = plotter, aes( x = rho, y= pi, col = source))+
  geom_line()+
  ylab(expression(pi))+
  xlab('Distance to Exon (Rho)')+
  scale_color_manual('', values = cbPalette)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17, angle = 0, vjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13)
  )


