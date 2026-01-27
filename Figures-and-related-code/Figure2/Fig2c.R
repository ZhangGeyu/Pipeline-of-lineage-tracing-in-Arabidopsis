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

library(phangorn)
library(dplyr)

# input the cell lineage tree generate in Fig2a
tree <- read.tree('ConsensusSeq_Tree_Plant1.nwk')
length(tree$tip.label)
traits <- data.frame(species = tree$tip.label)

# Remove germline sequences.
tree <- drop.tip(tree, tip = traits %>% filter(grepl("P", species)) %>% pull(species))
length(tree$tip.label)

# Remove rosette sequences.
Rosette_tip <- traits[grep("^R", traits$species), ]
tree <- drop.tip(tree, tip = Rosette_tip)
length(tree$tip.label)

tree <- makeNodeLabel(tree, method = "number", prefix = "N")
tree <- root(tree, outgroup = "reference", resolve.root = TRUE)
tree$edge.length[tree$edge.length == 0] <- 5e-10

# Classify the progeny and different branches.
traits$trait_value <- ifelse(grepl('P', traits$species), '0',
                             ifelse(substr(traits$species, 2, 2) == '2', '2',
                                    ifelse(substr(traits$species, 2, 2) == '3', '3',
                                           ifelse(substr(traits$species, 2, 2) == '4', '4',
                                                  ifelse(substr(traits$species, 2, 2) == '1', '1', '0')))))
head(traits)

rownames(traits) <- traits$species
traits$species <- NULL

trait_mat <- as.matrix(traits)
rownames(trait_mat) <- rownames(traits)
head(trait_mat)

trait_phyDat <- phyDat(trait_mat, type = "USER", levels = c("0","1",'2','3','4'))
obs_score <- parsimony(tree, trait_phyDat)

# Perform permutation test.
set.seed(123)
n_permute <- 1000
null_scores <- numeric(n_permute)

for (i in 1:n_permute) {
  permuted_states <- sample(trait_mat)  # Scramble trait labels
  names(permuted_states) <- rownames(trait_mat)
  perm_phyDat <- phyDat(permuted_states, type = "USER", levels = c("0","1",'2','3','4'))
  null_scores[i] <- parsimony(tree, perm_phyDat)
}

p_value <- mean(null_scores <= obs_score)
df_score <- data.frame(score = null_scores)

S_M_test <- ggplot(df_score, aes(x = score)) +
  geom_histogram(bins = 50, fill = "grey", color = "black") +
  geom_vline(xintercept = obs_score, color = "red", linetype = "dashed", size = 0.5) +
  labs(title = "Plant1 Parental",x = "Parsimony Score",y = "Count") +
  scale_x_continuous(breaks = seq(0, 1000, by = 20)) + 
  scale_y_continuous(breaks = seq(0, 1000, by = 50)) + 
  annotate("text", x = obs_score + 10, y = max(table(null_scores)) + 100,
           label = paste0("Observed = ", obs_score,'\n','P value = ',p_value),
           color = "red", hjust = 0, size = 3) + 
  Figure_Theme
S_M_test
