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

library(dplyr)
library(stringr)

# Readouts used in Fig1d
Sample_UMI_list <- c("B2-CL1_494","RL1_929","RL1_1155","B1-1-CL1_200","RL4_140","RL4_98","RL4_57")
Hotspot_Plant1 <- c("1092_T_-1C","1114_A_G","1226_C_T","1335_C_T","445_C_T","513_C_G","53_C_-2GA","77_G_A","795_G_T","965_G_A")

# Raw file
Parental_mutation <- read.table("Parental_CallSNP_Plant1.txt", header = T)
Parental_mutation <- Parental_mutation[!Parental_mutation$mut_info %in% Hotspot_Plant1,] # remove hotspot mutations

# Select the used Readouts from the dataframe
UMI_mut_example <- Parental_mutation %>%
  filter(SampleName_UMI %in% Sample_UMI_list) %>%
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

UMI_mut_example$y <- 0

UMI_mut_example$SampleName_UMI <- factor(UMI_mut_example$SampleName_UMI, levels = Sample_UMI_list)
UMI_mut_example$type <- factor(UMI_mut_example$type, levels = c('C>T/G>A','C>G/G>C','C>A/G>T','indel'))

# Plant 1 was analyzed using a different reference sequence than plants 2 and 3.
# Mutation positions were aligned to correspond with the mutation positions in plants 2 and 3.
UMI_mut_example$pos = UMI_mut_example$pos - 140
df_line <- data.frame(x = 1:(1353-140), y = rep(0, (1353-140)))

mut_Example <- ggplot() +
  geom_line(data = df_line, aes(x = x, y = y)) +
  geom_point(data = UMI_mut_example, aes(x = pos, y = y, colour = type)) +
  scale_colour_manual(values = c('#e63946','#30343f','#00a8e8','#4bad32')) + 
  scale_x_continuous(breaks = seq(0, 1300, by =100)) + 
  scale_y_continuous(breaks = seq(0, 0, by =0)) + 
  facet_wrap(~SampleName_UMI, ncol = 1, strip.position = "left") + 
  Figure_Theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y.left = element_text(angle = 0, hjust = 0),
        strip.background = element_blank()) 
mut_Example


