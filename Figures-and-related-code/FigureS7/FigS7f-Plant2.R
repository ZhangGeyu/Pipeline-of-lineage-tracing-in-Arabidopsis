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

HighFreqMut_list_886copy <- c('1213_C_+5GTGTG','886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC')

Parental_SNP <- read.table('Parental_CallSNP_Plant2_886copy.txt', header = T)
Parental_SNP <- Parental_SNP[!Parental_SNP$mut_info %in% HighFreqMut_list_886copy,]   # remove copy mutations
head(Parental_SNP)

UMI_mut_counts <- table(Parental_SNP$SampleName_UMI)
UMI_mut_counts <- as.data.frame(UMI_mut_counts)
names(UMI_mut_counts) <- c("SampleName_UMI", "count")
head(UMI_mut_counts)

Mut_count_distribution <- ggplot(UMI_mut_counts, aes(x = count)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 18,fill = "lightblue", colour = 'black', alpha = 0.5) +
  labs(x = "Mutation Count per Readout",y = "Density",
       title = "Plant2: mutation num in readouts") +
  scale_x_continuous(breaks = seq(0, 200, by = 2)) + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.1)) + 
  Figure_Theme
