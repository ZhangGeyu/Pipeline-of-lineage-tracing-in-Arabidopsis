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


# Raw file: 
Parental_mutation <- read.table("ConsensusSequence_CallSNP_Raw.txt", header = T)
head(Parental_mutation)

Parental_mutation <- Parental_mutation %>%
  mutate(
    # Define mutation types
    type = str_to_upper(paste0(ref, ">", alt))
  ) %>%
  mutate(
    type = case_when(
      str_detect(alt, "\\+|\\-") ~ "indel",
      type %in% c("C>T", "G>A") ~ "C>T/G>A",
      type %in% c("C>G", "G>C") ~ "C>G/G>C",
      type %in% c("C>A", "G>T") ~ "C>A/G>T",
      type %in% c("T>C", "A>G") ~ "T>C/A>G",
      type %in% c("T>A", "A>T") ~ "T>A/A>T",
      type %in% c("T>G", "A>C") ~ "T>G/A>C",
      TRUE ~ type
    )
  )

Mut_Type_Count <- as.data.frame(table(Parental_mutation$type))
colnames(Mut_Type_Count) <- c("type", "count")
Mut_Type_Count <- Mut_Type_Count[Mut_Type_Count$type != 'indel',]
Mut_Type_Count$type <- factor(Mut_Type_Count$type, levels = c('C>T/G>A','C>G/G>C','C>A/G>T',
                                                              'T>C/A>G','T>A/A>T','T>G/A>C'))

Mut_Type_Count_bar <- ggplot(Mut_Type_Count, 
                             aes(x = type, y = count/10000000)) + 
  geom_bar(stat = 'identity', width = 0.8, fill = 'grey', colour = 'black') + 
  xlab('Mutation type') + ylab('# of Mutations (x10^7)') + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 2, by = 0.25)) + 
  Figure_Theme + 
  theme(axis.text.x = element_text(angle=60, hjust=1))
