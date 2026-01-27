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

library(UpSetR)
library(grid)

SampleName <- 'B2-6-P1-1'   # Specify sample name to analyze
Sample_Haplotype <- read.table(paste0(SampleName, ".txt"), header = TRUE)
head(Sample_Haplotype)

n_labels <- length(unique(Offspring_CopyMut_sample$mut_info))
    
Progeny_Mut_list <- split(
  as.character(Sample_Haplotype$SampleName_UMI),
  as.character(Sample_Haplotype$mut_info))
  
# Use UpSetR to create Upset plot of mutation combinations
Progeny_Mut_Upset <- upset(
  fromList(Progeny_Mut_list),
  nsets = length(Progeny_Mut_list),
  nintersects = 50,
  mb.ratio = c(0.5, 0.5),
  order.by = "freq")
Progeny_Mut_Upset
