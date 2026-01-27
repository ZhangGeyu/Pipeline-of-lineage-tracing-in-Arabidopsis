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

library(dplyr)
library(stringr)

# Count mutations per UMI sample
UMI_mut_counts <- table(Parental_SNP$SampleName_UMI)
UMI_mut_counts <- as.data.frame(UMI_mut_counts)
names(UMI_mut_counts) <- c("SampleName_UMI", "count")
UMI_mut_counts <- UMI_mut_counts %>%
  mutate(SampleName = str_extract(as.character(SampleName_UMI), "^[^_]+"))
UMI_mut_counts$Group <- substr(UMI_mut_counts$SampleName,1,2)

head(UMI_mut_counts)
# Convert sample name to factor and specify level order
UMI_mut_counts$SampleName <- factor(UMI_mut_counts$SampleName, levels = c("B1-CL1","B1-1-CL1","B1-1-CL2","B1-1-CL3","B1-CL2","B1-2-CL1","B1-2-CL2","B1-2-CL3","B1-3-CL1","B1-3-CL3",
                                                                          "B1-4-CL2","B2-CL1","B2-CL3","B2-3-CL1","B2-3-CL3","B2-4-CL3","B2-CL6","B2-6-CL1","B3-5-CL2","B3-CL2",
                                                                          "B3-CL6","RL8","RL9","RL11","RL12","RL13"))
library(forcats)
Sample_UMI_MutCount_violin <- ggplot(UMI_mut_counts, 
                                     aes(x = fct_rev(SampleName), y = log10(count))) + 
  geom_violin(aes(fill = Group, colour = Group)) + 
  coord_flip() + 
  scale_fill_manual(values = c('#1ec0ff','#9381ff','#3ac569','#585858')) + 
  scale_colour_manual(values = c('#1ec0ff','#9381ff','#3ac569','#585858')) + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) + 
  xlab('SampleName') + ylab('log10(# of mutations per readout)') + 
  ggtitle('Plant3') + 
  Figure_Theme
