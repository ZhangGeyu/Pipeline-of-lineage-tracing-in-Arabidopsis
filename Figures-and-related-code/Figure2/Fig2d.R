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

Correlation_group <- read.table('Sample_MutFreq_correlation.txt',header = TRUE)
Correlation_group$group <- factor(Correlation_group$group,levels = c('Within_Branch','Between_Branch'))
head(Correlation_group)

df_counts <- Correlation_group %>% group_by(group) %>% summarise(n = n(), y_pos = 0.98) 

Correlation_violin <- ggplot(Correlation_group, 
                             aes(x = group, y = pcc)) +
  geom_boxplot(width = 0.6, colour = '#333533',fill = 'grey',
               notch = TRUE, outlier.shape = NA) +
  geom_jitter(size = 0.5, shape = 16, width = 0.3) +
  geom_signif(comparisons = list(c('Within_Branch','Between_Branch')), 
              color = "black", textsize = 2.5,
              tip_length = 0.005,y_position = 1, test = 'wilcox.test') +
  ylab("Pearson's Correlation Coefficient") + 
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.1)) + 
  geom_text(data = df_counts, aes(x = group, y = y_pos, label = paste0("n=", n)), 
            color = 'black', size = 3) +
  Figure_Theme
Correlation_violin
