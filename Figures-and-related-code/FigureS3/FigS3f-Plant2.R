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

# Plant 2

MutFreq_In_ProgenySample <- read.table('MutFreq_In_ProgenySample_Plant2_886copy.txt', header = TRUE, sep = "\t")

MutFreq_In_ProgenySample <- MutFreq_In_ProgenySample %>% filter(mut_freq >= 0.5)  # keep records with mutation frequency ≥0.5
Progeny_MutCount <- MutFreq_In_ProgenySample %>% group_by(SampleName) %>% summarise(count = n()) # Group by sample and count number of mutations

Progeny_MutCount$SampleName <- factor(Progeny_MutCount$SampleName,
                                      levels = c("B1-P2-2","B1-P2-8","B1-P2-11","B1-P2-13","B1-P2-20","B1-P3-10","B1-P3-11","B1-P3-13","B1-P3-14",
                                                 "B1-1-P1-4","B1-1-P1-7","B1-1-P1-9","B1-1-P1-11","B1-1-P1-13","B1-1-P1-14","B1-1-P2-4","B1-1-P2-6","B1-1-P2-7","B1-1-P2-9","B1-1-P2-11","B1-1-P2-12","B1-1-P2-15","B1-1-P2-18","B1-2-P2-9","B1-2-P2-12","B1-2-P3-4","B1-2-P3-9","B1-2-P3-12","B1-2-P3-13",
                                                 "B2-P2-2","B2-P2-3","B2-P2-9","B2-P3-1","B2-P3-3","B2-P3-5","B2-P3-8","B2-P3-9",
                                                 "B2-1-P1-10","B2-1-P1-11","B2-1-P2-8","B2-2-P1-4","B2-2-P1-6","B2-2-P1-12","B2-2-P2-6","B2-2-P2-9",
                                                 "B3-P2-3","B3-P2-6","B3-P3-2","B3-P3-9","B3-P3-11",
                                                 "B4-P2-1","B4-P2-3","B4-P2-6","B4-P2-10","B4-P2-13"))

Progeny_MutCount_bar <- ggplot(Progeny_MutCount,aes(x = fct_rev(SampleName), y = count)) + 
  geom_bar(stat = 'identity', colour = 'black', fill = 'grey') +
  coord_flip() + 
  scale_y_continuous(breaks = seq(0, 500, by = 5)) + 
  xlab('SampleName') + ylab('# of mutations per haplotype') + 
  ggtitle('Plant2') + 
  Figure_Theme
