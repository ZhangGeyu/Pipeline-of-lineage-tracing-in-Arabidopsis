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

# Plant3

MutFreq_In_ParentalSample <-  read.table('MutFreq_In_ParentalSample_Plant3_181copy.txt', header=T)
MutFreq_In_ParentalSample$mut_freq <- as.numeric(MutFreq_In_ParentalSample$mut_freq)
MutFreq_In_ParentalSample$branch <- substr(MutFreq_In_ParentalSample$SampleName, 2, 2)
head(MutFreq_In_ParentalSample)

Plant_UMI_num <- sum(unique(MutFreq_In_ParentalSample[, c("SampleName", "Sample_count")])$Sample_count)
Plant_UMI_num

# Calculate the mutation frequency in the whole plant

MutFreq_In_Plant <- MutFreq_In_ParentalSample %>% 
  # Count the total number of UMIs with each mutation.
  group_by(mut_info) %>% 
  summarise(total_mut_count = sum(mut_info_count),
            total_reads = sum(Sample_count), 
            sample_count = n(),.groups = 'drop') %>%
  mutate(mut_freq = total_mut_count /Plant_UMI_num) %>%
  select(mut_info,total_mut_count,sample_count,mut_freq)
head(MutFreq_In_Plant)


# extract branch information for each mutation type from the raw data.
branch_info <- MutFreq_In_ParentalSample %>%
  distinct(mut_info, branch) %>%  # Get unique mutation-branch combinations
  group_by(mut_info) %>%          # Group by mutation types
  summarise(branch_count = n())   # Count the number of different branches
head(branch_info)

# Merge the branch information into the result data frame
MutFreq_In_Plant <- MutFreq_In_Plant %>%
  left_join(branch_info, by = "mut_info") %>%
  mutate(branch_count = ifelse(is.na(branch_count), 0, branch_count))
head(MutFreq_In_Plant)


# Create a new column: log10 mutation frequency (handling zero values issue).
MutFreq_In_Plant$log10_mut_freq <- log10(MutFreq_In_Plant$mut_freq)   # Add a small value to avoid log(0).
MutFreq_In_Plant<-subset(MutFreq_In_Plant,MutFreq_In_Plant$mut_freq!=1)

# Create bins for log10_mut_freq (width 0.1)
MutFreq_In_Plant$log10_mut_freq_bin <- cut(MutFreq_In_Plant$log10_mut_freq, 
                                 breaks = seq(floor(min(MutFreq_In_Plant$log10_mut_freq, na.rm = TRUE)), 
                                              ceiling(max(MutFreq_In_Plant$log10_mut_freq, na.rm = TRUE)), 
                                              by = 0.1),
                                 include.lowest = TRUE)
head(MutFreq_In_Plant)

# Calculate outliers for sample_count within each bin
MutFreq_In_Plant <- MutFreq_In_Plant %>%
  group_by(log10_mut_freq_bin) %>%
  mutate(
    Q1 = quantile(sample_count, 0.25, na.rm = TRUE),
    Q3 = quantile(sample_count, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    Down_outlier = sample_count < lower_bound,
    Up_outlier = sample_count > upper_bound
  ) %>%
  ungroup()

head(MutFreq_In_Plant)

Hotspot_mut_identify <- ggplot(MutFreq_In_Plant, aes(x = log10_mut_freq_bin, y = sample_count)) +
  geom_boxplot(outlier.shape = NA,alpha = 0.3, color = "black") +
  geom_point(data = filter(MutFreq_In_Plant, Up_outlier & branch_count >= 4),
             position = position_jitter(width = 0.2, height = 0.2),
             aes(color = factor(branch_count)),size = 1.5,alpha = 0.8,shape = 16) +
  geom_point(data = filter(MutFreq_In_Plant, !Up_outlier),
             position = position_jitter(width = 0.2, height = 0.2),
             color = 'darkgrey',size = 0.5, alpha = 0.5,shape = 16) +
  scale_color_manual(values = c("1" = "#52b788","2" = "#1e88e5",
                                "3" = "#7b2cbf","4" = "#e63946", "5" = "#e63946")) +
  scale_y_continuous(breaks = seq(0, 50, by = 5)) + 
  labs(x = "Log10(Mutation Frequency) (binned by 0.1)", 
       y = "Number of Samples with Mutation",
       title = "Plant3") +
  Figure_Theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hotspot_mutation <- unique(subset(MutFreq_In_Plant, Up_outlier & branch_count >= 4)$mut_info)
hotspot_mutation
