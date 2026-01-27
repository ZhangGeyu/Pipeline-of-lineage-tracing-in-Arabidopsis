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

## Clustered Progeny

library(phangorn)
# Read phylogenetic tree file
tree <- read.tree('ConsensusSeq_Tree_Plant2.nwk')
traits <- data.frame(species = tree$tip.label)
# Label internal nodes of the tree
tree <- makeNodeLabel(tree, method = "number", prefix = "N")
# Root the tree using the reference sequence as outgroup
tree <- root(tree, outgroup = "reference", resolve.root = TRUE)
# Replace branches with length 0 with a tiny positive number
tree$edge.length[tree$edge.length == 0] <- 5e-10

Offspring_1 <- c('B1_P3_10_progeny','B1_2_P2_9_progeny','B1_1_P2_4_progeny','B1_1_P1_9_progeny',
                 'B2_2_P1_6_progeny','B1_1_P2_12_progeny','B2_P2_9_progeny','B3_P3_11_progeny',
                 'B2_2_P1_12_progeny','B4_P2_1_progeny','B2_2_P2_6_progeny')
traits$trait_value <- ifelse(traits$species %in% Offspring_1, 1, 0)

rownames(traits) <- traits$species
traits$species <- NULL
# Convert traits to matrix
trait_mat <- as.matrix(traits)
rownames(trait_mat) <- rownames(traits)
head(trait_mat)
# Convert trait matrix to phyDat format
trait_phyDat <- phyDat(trait_mat, type = "USER", levels = c("0","1"))
obs_score <- parsimony(tree, trait_phyDat)
# Perform permutation test
set.seed(123)
n_permute <- 1000
null_scores <- numeric(n_permute)

for (i in 1:n_permute) {
  # Randomly permute trait states
  permuted_states <- sample(trait_mat) 
  names(permuted_states) <- rownames(trait_mat)
  perm_phyDat <- phyDat(permuted_states, type = "USER", levels = c("0","1"))
  # Calculate parsimony score after permutation
  null_scores[i] <- parsimony(tree, perm_phyDat)
}

p_value <- mean(null_scores <= obs_score)
df_score <- data.frame(score = null_scores)

S_M_test <- ggplot(df_score, aes(x = score)) +
  geom_histogram(bins = 10, fill = "grey", color = "black") +
  geom_vline(xintercept = obs_score, color = "red", linetype = "dashed", size = 0.5) +
  labs(title = "Plant2 Cluster Progeny",x = "Parsimony Score",y = "Count") +
  scale_x_continuous(breaks = seq(0, 1000, by = 1)) + 
  scale_y_continuous(breaks = seq(0, 1000, by = 200)) + 
  annotate("text", x = obs_score + 0.5, y = max(table(null_scores)) - 100,
           label = paste0("Observed = ", obs_score,'\n','P value = ',p_value),
           color = "red", hjust = 0, size = 3) + 
  Figure_Theme
S_M_test
