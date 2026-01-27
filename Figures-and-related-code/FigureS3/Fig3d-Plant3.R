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

# Plant3

Hotspot_181copy <- c("1007_C_T","1013_C_T","1100_G_A","1186_A_G","1242_G_A","306_T_C","77_C_T","850_C_T","941_C_T")

MutFreq_In_ProgenySample <- read.table('MutFreq_In_ProgenySample_Plant3_181copy.txt', header = T)

MutFreq_In_ProgenySample <- MutFreq_In_ProgenySample[(MutFreq_In_ProgenySample$mut_freq >= 0.05)&(!MutFreq_In_ProgenySample$mut_info %in% Hotspot_181copy),]
head(MutFreq_In_ProgenySample)

Mut_freq_distribution <- ggplot(MutFreq_In_ProgenySample, aes(x = mut_freq)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(x = "Mutation frequency",y = "Density",
       title = "Plant3: mutation frequency in progeny samples") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) + 
  Figure_Theme
Mut_freq_distribution
