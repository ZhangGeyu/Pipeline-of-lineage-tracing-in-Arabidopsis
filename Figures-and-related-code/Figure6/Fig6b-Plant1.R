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

PassToOffspring_df <- read.table(paste0(path,'Mut_PassToOffspring_BranchShared-Plant1.txt'), header = TRUE)
head(PassToOffspring_df)

Parent_MutCount_bar <- ggplot(PassToOffspring_df,
                              aes(x = BranchCount, y = Parent_Mut, fill = BranchCount)) + 
  geom_bar(stat = 'identity',colour = 'black') + 
  scale_y_continuous(breaks = seq(0, 10000, by = 500)) + xlab('') + 
  scale_fill_manual(values = c('#ade8f4','#48cae4','#0096c7','#0077b6')) + 
  Figure_Theme +
  theme(axis.text.x = element_blank(),
        legend.position="none")
Parent_MutCount_bar

PassToOffspring_MutCount_bar <- ggplot(PassToOffspring_df,
                                       aes(x = BranchCount, y = PassToOffspring, fill = BranchCount)) + 
  geom_bar(stat = 'identity',colour = 'black') + 
  scale_y_continuous(breaks = seq(0, 500, by = 50)) + xlab('') + 
  scale_fill_manual(values = c('#ade8f4','#48cae4','#0096c7','#0077b6')) + 
  Figure_Theme +
  theme(axis.text.x = element_blank(),
        legend.position="none")
PassToOffspring_MutCount_bar

PassToOffspring_MutFraction_bar <- ggplot(PassToOffspring_df,
                                          aes(x = BranchCount, y = Fraction, fill = BranchCount)) + 
  geom_bar(stat = 'identity',colour = 'black') + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) + xlab('') + 
  scale_fill_manual(values = c('#ade8f4','#48cae4','#0096c7','#0077b6')) + 
  Figure_Theme +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
        legend.position="none")
PassToOffspring_MutFraction_bar
