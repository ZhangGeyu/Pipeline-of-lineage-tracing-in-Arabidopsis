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

Bootstrap <- read.table('Tree_bootstrap.txt', header = TRUE)

Bootstrap_sub <- subset(Bootstrap, bootstrap != 0)
median_val <- median(Bootstrap_sub$bootstrap, na.rm = TRUE)

Bootstrap_Hist <- ggplot(Bootstrap_sub, aes(x = bootstrap)) +
  geom_histogram(bins = 50, width = 1, colour = 'black', fill = 'grey') +
  geom_vline(xintercept = median_val,
             colour = "red", linetype = "dashed", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  scale_y_continuous(breaks = seq(0, 500, by = 10)) + 
  annotate("text", x = median_val -0.5, y = 50,
           label = paste0("Median = ", median_val),
           color = "red", hjust = 0, size = 3) + 
  Figure_Theme
