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

distance <- read.delim2('different_mutation_between_sample_final.txt',stringsAsFactors=F,header=F)
sorted_pairs <- apply(distance[c("V1", "V2")], 1, function(x) paste(sort(x), collapse =":"))
distance$V3<- as.numeric(distance$V3)
distance$V4<- as.numeric(distance$V4)

model <- lm(V3 ~ V4, data = distance)
intercept <- coef(model)[1]  
slope <- coef(model)[2]     
r_squared <- summary(model)$r.squared
eq_label <- sprintf("y = %.2f + %.2fx", intercept, slope)
r2_label <- sprintf("R² = %.2f", r_squared)

ggplot(distance, aes(x = V4, y = V3)) +
  geom_point(color = "#2c7fb8", size = 3, alpha = 0.7) + 
  geom_smooth(method = "lm", formula = y ~ x,
              color = "#e41a1c", fill = "#fed976",
              se = TRUE, level = 0.95) +
  labs(x = "Number of internodes between two cauline leaves", 
       y = "Number of unshared mutations between two cauline leaves") +
  stat_cor(label.x = 6, label.y =160) +
  annotate("text", x=6, y = 150,
           label = paste(eq_label),
           size = 4, hjust = 0)+Figure_Theme
