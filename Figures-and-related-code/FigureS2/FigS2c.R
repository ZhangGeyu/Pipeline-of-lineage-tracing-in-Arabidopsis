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

Allele_Redundancy_Count <- read.table('Allele_Redundancy_Count.txt', sep = '\t', header = TRUE)
Allele_Redundancy_Count$Redundancy <- factor(Allele_Redundancy_Count$Redundancy,levels = c('1','2','3','4','>=5'))
head(Allele_Redundancy_Count)

Allele_Redundancy_Count_bar <- ggplot(Allele_Redundancy_Count, 
                                      aes(x = Redundancy, y = count/1000)) + 
  geom_bar(stat = 'identity', width = 0.8, colour = 'black', fill = 'grey') + 
  xlab('Redundancy level') + ylab('count (×10^3)') + 
  scale_y_continuous(breaks = seq(0, 5, by = 1)) + 
  Figure_Theme
