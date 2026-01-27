Figure_Theme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  theme(plot.title=element_text(size=8))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5))+
  theme(axis.text.x = element_text(colour="black", size=8))+
  theme(axis.text.y = element_text(colour="black", size=8))+
  theme(axis.title=element_text(size=8))+
  theme(axis.ticks=element_line(colour="black",size=0.5))+
  theme(axis.line=element_line(colour="black")) + 
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=8))+
  theme(axis.line=element_blank()) + 
  theme(panel.border = element_rect(fill=NA, size=1))

DistanceFromTargetSite <- read.table('DistanceFromTargetSite.txt', sep = '\t', header = TRUE)
head(DistanceFromTargetSite)

Distance_Density <- ggplot(DistanceFromTargetSite, 
                           aes(x = Distance)) + 
  scale_x_continuous(limits = c(-20,40),breaks = seq(-20, 40, by = 20)) + 
  geom_density(adjust = 1) + 
  xlab('Distance to binding motif (bp)') + 
  Figure_Theme
Distance_Density
