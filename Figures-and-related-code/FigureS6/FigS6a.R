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

# somatic mutation spectrum

library(tidyverse)

Self1 <- c('#E72726', '#060709', '#1EBEF0', '#A1CF63', '#CACACA', '#EEC8C5')

read_and_process_data <- function(file_path, sample_name) {
  data <- read.delim2(file_path, stringsAsFactors = FALSE, header = FALSE)
  colnames(data) <- c('chr', 'pos', 'mut', 'mut_number', 'total', 'fre')

  data <- data %>%
    mutate(
      total = as.numeric(total),
      fre = as.numeric(fre),
      mut_number = as.numeric(mut_number),
      sample = sample_name
    ) %>%
    filter(fre != 1,
           mut_number >= 3,
           total >= 30,
           total <= 100,
           !chr %in% c('ChrC', 'ChrM', 'XF675_1kb', 'XF4363_AI'))
  
  return(data)
}

file_paths <- c(
  WT = '/../UMIC_vcf_test_WT.txt',    
  plant1 = '/../UMIC_vcf_test_plant1.txt',            
  plant2 = '/../UMIC_vcf_test_plant2.txt',   
  plant3 = '/../UMIC_vcf_test_plant3.txt'  
)

all_data <- map2_dfr(file_paths, names(file_paths), read_and_process_data)

# Filter out mutations that appear in only a single sample
unique_mutations <- all_data %>%
  add_count(chr, pos, mut) %>%
  filter(n == 1) %>%
  select(-n)

classify_mutation_spectrum <- function(data) {
  data %>%
    mutate(
      spectrum = case_when(
        mut %in% c("C>T", "G>A") ~ "C>T/G>A",
        mut %in% c("C>G", "G>C") ~ "C>G/G>C",
        mut %in% c("C>A", "G>T") ~ "C>A/G>T",
        mut %in% c("T>C", "A>G") ~ "T>C/A>G",
        mut %in% c("T>A", "A>T") ~ "T>A/A>T",
        mut %in% c("T>G", "A>C") ~ "T>G/A>C",
        TRUE ~ NA_character_
      )
    )
}

create_spectrum_plot <- function(data, sample_name) {
  sample_data <- data %>%
    filter(sample == sample_name) %>%
    classify_mutation_spectrum()
  
  spectrum_counts <- sample_data %>%
    count(spectrum, name = "Freq") %>%
    filter(!is.na(spectrum)) %>%
    mutate(spectrum = factor(spectrum, 
                             levels = c("C>T/G>A", "C>G/G>C", "C>A/G>T",
                                        "T>C/A>G", "T>A/A>T", "T>G/A>C")))

  all_levels <- c("C>T/G>A", "C>G/G>C", "C>A/G>T",
                  "T>C/A>G", "T>A/A>T", "T>G/A>C")
  spectrum_counts <- spectrum_counts %>%
    complete(spectrum = factor(all_levels, levels = all_levels), 
             fill = list(Freq = 0))

  ggplot(spectrum_counts, aes(x = spectrum, y = Freq, fill = spectrum)) +
    geom_bar(stat = 'identity', color = NA, width = 0.8) +
    scale_fill_manual(values = Self1) +
    labs(y = "Number of mutations", x = '') +
    ggtitle(paste("Mutation Spectrum -", sample_name)) +
    Figure_Theme +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
}

samples <- c("WT", "plant1", "plant2", "plant3")
plots <- map(samples, ~ create_spectrum_plot(unique_mutations, .x))
names(plots) <- samples

walk2(plots, names(plots), ~ {
  print(.x)
  cat("plot:", .y, "\n")
})

