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

# plant2

HighFreqMut_list_886copy <- c('1213_C_+5GTGTG','886_G_A','862_G_C','874_G_A','837_G_C','904_G_T','841_C_T','458_T_-20AACAGGGTAATGAGCCGCAC')
Hotspot_886copy <- c("")

Parental_SNP <- read.table('Parental_CallSNP_Plant2_886copy.txt', header = T)
Parental_SNP <- Parental_SNP[(!Parental_SNP$mut_info %in% HighFreqMut_list_886copy)&(!Parental_SNP$mut_info %in% Hotspot_886copy),]        # Remove copy mutations and hotspots
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

unique(UMI_mut_counts$SampleName)
# Convert sample name to factor and specify level order
UMI_mut_counts$SampleName <- factor(UMI_mut_counts$SampleName, levels = c("B1-CL1","B1-1-CL1","B1-1-CL2","B1-2-CL1","B1-3-CL1","B1-CL2","B1-CL3","B2-CL1","B2-1-CL1","B2-1-CL2",
                                                                          "B2-CL2","B2-2-CL1","B2-2-CL2","B2-CL3","B2-3-CL1","B2-3-CL2","B3-CL1","B3-1-CL1","B3-1-CL2","B3-CL2",
                                                                          "B3-2-CL1","B3-CL3","B3-3-CL1","B3-3-CL2","B4-CL1","B4-1-CL1","B4-CL2","B4-2-CL1","B4-2-CL2","RL1",
                                                                          "RL2","RL3","RL4","RL6","RL7","RL8","RL9"))
library(forcats)
Sample_UMI_MutCount_violin <- ggplot(UMI_mut_counts, 
                                     aes(x = fct_rev(SampleName), y = log10(count))) + 
  geom_violin(aes(fill = Group, colour = Group)) + 
  coord_flip() + 
  scale_fill_manual(values = c('#1ec0ff','#9381ff','#3ac569','#ee9b00','#585858')) + 
  scale_colour_manual(values = c('#1ec0ff','#9381ff','#3ac569','#ee9b00','#585858')) + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) + 
  xlab('SampleName') + ylab('log10(# of mutations per readout)') + 
  ggtitle('Plant2') + 
  Figure_Theme
