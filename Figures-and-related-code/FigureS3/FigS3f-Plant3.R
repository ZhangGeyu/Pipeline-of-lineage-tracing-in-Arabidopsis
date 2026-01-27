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

Haplotype_ProgenySample <- read.table('Offspring_haplotype.txt', header = TRUE, sep = "\t")
Progeny_MutCount <- Haplotype_ProgenySample %>% group_by(SampleName_UMI) %>% summarise(count = n()) # Group by sample and count number of mutations
head(Progeny_MutCount)

Progeny_MutCount$SampleName_UMI <- factor(Progeny_MutCount$SampleName_UMI,
                                      levels = c("B1-P1-1_1st","B1-P1-2_1st","B1-P1-3_1st","B1-P1-3_2nd","B1-P2-1_1st","B1-P2-2_1st","B1-P2-3_1st","B1-P2-4_1st","B1-P2-5_1st","B1-P6-1_1st",
                                                 "B1-P6-2_1st","B1-P6-3_1st","B1-P6-3_2nd","B1-P6-5_1st","B1-4-P1-2_1st","B1-4-P1-3_1st","B1-4-P1-3_2nd","B1-4-P1-4_1st","B1-4-P2-1_1st","B1-4-P2-1_2nd",
                                                 "B1-4-P2-5_1st","B1-4-P4-2_1st","B1-4-P4-3_1st","B1-4-P4-3_2nd","B2-P5-1_1st","B2-P5-3_1st","B2-6-P1-1_1st","B2-6-P1-1_2nd","B2-6-P1-2_1st","B2-6-P1-3_1st",
                                                 "B2-6-P1-6_1st","B2-6-P2-1_1st","B2-6-P2-1_2nd","B2-6-P2-3_1st","B2-6-P2-4_1st","B2-6-P2-5_1st","B2-6-P5-1_1st","B2-6-P5-2_1st","B2-6-P5-4_1st","B2-6-P5-4_2nd",
                                                 "B3-P2-1_1st","B3-P3-2_1st","B3-P6-1_1st","B3-P6-2_1st","B3-P6-4_1st","B3-P6-5_1st"))

Progeny_MutCount_bar <- ggplot(Progeny_MutCount,aes(x = fct_rev(SampleName_UMI), y = count)) + 
  geom_bar(stat = 'identity', colour = 'black', fill = 'grey') +
  coord_flip() + 
  scale_y_continuous(breaks = seq(0, 500, by = 5)) + 
  xlab('SampleName') + ylab('# of mutations per haplotype') + 
  ggtitle('Plant3') + 
  Figure_Theme
