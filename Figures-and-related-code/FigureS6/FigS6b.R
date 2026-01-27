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

genome<-  read.delim2('genome_coverage.txt',stringsAsFactors=F,header=F,sep =' ')
ggplot(genome, aes(x = V2, y = V3)) +
  geom_line(color = "black", size = 0.8) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) +
  labs(title = "Genomic region",
       x = "Coordinate on Chromosome 1 reference (bp)",
       y = "Number of reads") +
  Figure_Theme + 
  theme_minimal()

AI<-  read.delim2('AI_coverage.txt',stringsAsFactors=F,header=F,sep =' ') 
ggplot(AI, aes(x = V2, y = V3)) +
  geom_line(color = "black", size = 0.8) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) +
  labs(title = "AI region",
       x = "Coordinate of the AI plasmid reference (bp)",
       y = "Number of reads") +
  Figure_Theme + 
  theme_minimal()

readout <-  read.delim2('EGFP-1kb_depth.txt',stringsAsFactors=F,header=F,sep =' ')
ggplot(readout, aes(x = V2, y = V3)) +
  geom_line(color = "black", size = 0.8) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) +
  labs(title = "Readout region",
       x = "Coordinate of the readout plasmid reference (bp)",
       y = "Number of reads") +
  Figure_Theme + 
  theme_minimal()
