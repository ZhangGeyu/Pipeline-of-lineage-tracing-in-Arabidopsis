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

# Plant2

Hotspot_886copy <- c("")
MutFreq_In_ProgenySample <- read.table('MutFreq_In_ProgenySample_Plant2_886copy.txt', header = T)

MutFreq_In_ProgenySample <- MutFreq_In_ProgenySample[(MutFreq_In_ProgenySample$mut_freq >= 0.1)&(!MutFreq_In_ProgenySample$mut_info %in% Hotspot_886copy),]
head(MutFreq_In_ProgenySample)

Mut_freq_distribution <- ggplot(MutFreq_In_ProgenySample, aes(x = mut_freq)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(x = "Mutation frequency",y = "Density",
       title = "Plant2: mutation frequency in progeny samples") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) + 
  Figure_Theme
Mut_freq_distribution
