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

# Raw file
Parental_SNP <- read.table('Parental_CallSNP_Plant1.txt', header = T)

Parental_SNP <- Parental_SNP[!Parental_SNP$mut_info %in% Hotspot_Plant1,]   # Remove hotspot mutations
head(Parental_SNP)

library(dplyr)
library(stringr)

UMI_mut_counts <- table(Parental_SNP$SampleName_UMI)
UMI_mut_counts <- as.data.frame(UMI_mut_counts)
names(UMI_mut_counts) <- c("SampleName_UMI", "count")
UMI_mut_counts <- UMI_mut_counts %>%
  mutate(SampleName = str_extract(as.character(SampleName_UMI), "^[^_]+"))    #Extract sample name (portion before underscore) from UMI identifier
UMI_mut_counts$Group <- substr(UMI_mut_counts$SampleName,1,2)      # Create Group column based on branch

unique(UMI_mut_counts$SampleName)

UMI_mut_counts$SampleName <- factor(UMI_mut_counts$SampleName, levels = c("B1-CL1","B1-1-CL1","B1-1-CL2","B1-CL4","B1-CL5",
                                                                          "B2-CL1","B2-1-CL1","B2-CL2",
                                                                          "B3-CL1","B3-1-CL1","B3-CL2","B3-2-CL1","B3-2-CL2","B3-CL3","B3-CL4",
                                                                          "RL1","RL4","RL6","RL7","RL10","RL11"))

Sample_UMI_MutCount_violin <- ggplot(UMI_mut_counts, 
                                     aes(x = SampleName, y = log10(count))) + 
  geom_violin(aes(fill = Group, colour = Group)) + 
  scale_fill_manual(values = c('#1ec0ff','#9381ff','#3ac569','#585858')) + 
  scale_colour_manual(values = c('#1ec0ff','#9381ff','#3ac569','#585858')) + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) + 
  xlab('SampleName') + ylab('log10(# of mutations per readout)') + 
  ggtitle('Plant1') + 
  Figure_Theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
