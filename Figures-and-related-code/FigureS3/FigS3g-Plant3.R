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

Mut_Type_Count <- read.table('MutTypeCount_181copy_Post-germination.txt', sep = '\t', header = TRUE)
Mut_Type_Count$Mut_type <- factor(Mut_Type_Count$Mut_type, levels = c('C>T/G>A','C>G/G>C','C>A/G>T',
                                                                      'T>C/A>G','T>A/A>T','T>G/A>C'))
head(Mut_Type_Count)

Mut_Type_Count_bar <- ggplot(Mut_Type_Count, 
                             aes(x = Mut_type, y = count/10000)) + 
  geom_bar(stat = 'identity', width = 0.8) + 
  xlab('Mutation type') + ylab('# of Mutations (×10^4)') + 
  ggtitle('Plant3: Post-germination mutations') + 
  scale_y_continuous(breaks = seq(0, 3, by = 0.1)) + 
  Figure_Theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
