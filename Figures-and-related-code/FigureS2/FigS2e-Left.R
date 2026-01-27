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

Hotspot_Plant1 <- c("1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A")

Offspring_MutCount <- read.table('Offspring_MutCount_Plant1.txt', header = T)
Offspring_MutCount <- Offspring_MutCount[(Offspring_MutCount$mut_freq >= 0.5)&(!Offspring_MutCount$mut_info %in% Hotspot_Plant1),]    # Remove hotspot mutations

Offspring_MutCount <- Offspring_CallSNP %>% group_by(SampleName_UMI) %>% summarise(count = n())
head(Offspring_MutCount)

Offspring_MutCount$SampleName_UMI <- factor(Offspring_MutCount$SampleName_UMI,
                                            levels = c("B1-1-P4-2","B1-1-P9-2","B1-2-P1-2","B1-2-P6-1","B1-2-P8-4",
                                                       "B2-P1-1","B2-P1-2","B2-P2-1","B2-P2-2","B2-P4-1","B2-P5-1","B2-P5-2","B2-2-P1-1","B2-2-P3-1","B2-2-P3-2",
                                                       "B3-P3-2","B3-P3-3","B3-P3-4","B3-1-P2-1","B3-1-P2-4","B3-1-P2-6","B3-1-P4-1","B3-2-P4-1","B3-2-P5-2"))

Offspring_MutCount_hist <- ggplot(Offspring_MutCount,
                                  aes(x = SampleName_UMI, y = count)) + 
  geom_bar(stat = 'identity', colour = 'black', fill = 'grey') +
  scale_y_continuous(breaks = seq(0, 500, by = 20)) + 
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5))+
  theme(axis.line=element_blank()) + 
  theme(panel.border = element_rect(fill=NA, size=1))
Offspring_MutCount_hist
