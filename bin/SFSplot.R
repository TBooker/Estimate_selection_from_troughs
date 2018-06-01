rm(list=ls())

library(ggplot2)
library(grDevices)
library(scales)
library(reshape2)
cbPalette <- c("#56B4E9",  "#D55E00", "#E69F00")

Ns400 <- read.csv('~/PhD/Coding/Estimate_selection_from_troughs/simulations/NewSimulations/SFSplot/Ns400_pa0.000125.sfs.csv')
Ns5 <- read.csv('~/PhD/Coding/Estimate_selection_from_troughs/simulations/NewSimulations/SFSplot/Ns5_pa0.01.sfs.csv')

sfs <- rbind(Ns400,Ns5)


str(sfs)

sfsN <- sfs[sfs$Source != 'Synonymous',]

vnames <-list(
  'Ns400_pa0.000125'  =  expression(N[e]*'s = 400 ; '*p[a]*' =  0.000125'),
  'Ns5_pa0.01 ' =  expression(N[e]*'s = 5 ; '*p[a]*' =  0.01')
  )

vlabeller <- function(variable,value){
  return(vnames[value])
}

ggplot(data = sfs, aes(x=alleles, y = sfsCount, fill = Source))+
  geom_bar(stat = 'identity', position = 'dodge')+
  scale_y_sqrt('Number of Polymorphic Sites',labels = comma)+
  facet_grid(~label, labeller = vlabeller)+
#  facet_grid(Source~label)+ 
  scale_fill_manual('' ,values= cbPalette)+
  scale_x_continuous('Number of Derived Alleles', breaks = seq(1,19))+
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

