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

# plant3

HighFreqMut_list_181copy <- c('1213_C_+5GTGCT','181_G_A','182_G_A','183_G_A','184_G_A','736_G_A','601_G_A','657_G_A','719_C_T','238_G_A','658_G_A','663_C_T')
Hotspot_181copy <- c("1007_C_T","1013_C_T","1100_G_A","1186_A_G","1242_G_A","306_T_C","77_C_T","850_C_T","941_C_T")

Parental_SNP <- read.table('Parental_CallSNP_Plant3_181copy.txt', header = T)
Parental_SNP <- Parental_SNP[(!Parental_SNP$mut_info %in% HighFreqMut_list_181copy)&(!MutFreq_In_ProgenySample$mut_info %in% Hotspot_181copy),]        # Remove copy mutations and hotspots
head(Parental_SNP)

UMI_mut_counts <- table(Parental_SNP$SampleName_UMI)
UMI_mut_counts <- as.data.frame(UMI_mut_counts)
names(UMI_mut_counts) <- c("SampleName_UMI", "count")
head(UMI_mut_counts)

Mut_count_distribution <- ggplot(UMI_mut_counts, aes(x = count)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,fill = "lightblue", colour = 'black', alpha = 0.5) +
  labs(x = "Mutation Count per Readout",y = "Density",
       title = "Plant3: mutation num in readouts") +
  scale_x_continuous(breaks = seq(0, 200, by = 20)) + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.01)) + 
  Figure_Theme
