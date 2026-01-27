library(ggplot2)

Figure_Theme <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank()) + theme(plot.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 8)) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.title = element_text(size = 8)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.5)) +
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.title = element_text(size = 8)) +theme(legend.text = element_text(size = 8)) +
  theme(axis.line = element_blank()) +theme(panel.border = element_rect(fill = NA, size = 1))


create_plot <- function(data, plant_name, chi_text, show_axes = FALSE) {
  n_branches <- length(unique(data$branch))
  y_max <- max(data$number) * 1.15
  p <- ggplot(data, aes(x = branch, y = number, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("early" = "#beaed4", "late" = "#7fc97f")) +Figure_Theme
  p <- p + geom_text(aes(label = number),
                     position = position_dodge(width = 0.7),
                     vjust = -0.5,  size = 3)  
  
  p <- p + ylim(0, y_max)
  
  if (!show_axes) {
    p <- p + theme( legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), axis.line = element_blank())}
  p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))
  p <- p + annotate("text",x = n_branches,  y = y_max * 0.95, label = chi_text,
                    hjust = 1, vjust = 1, size = 2.5) 
    p <- p + annotate("text", x = 0.5,  y = y_max * 0.95, label = plant_name,
                    hjust = 0, vjust = 1,size = 6, fontface = "bold")
  return(p)}

calculate_chi_text <- function(values_matrix) {
  chi_test <- chisq.test(values_matrix)
  chi_value <- round(chi_test$statistic, 3)
  df_value <- chi_test$parameter
  p_value <- ifelse(chi_test$p.value < 0.001, "< 0.001",  round(chi_test$p.value, 3))
  paste0("χ² = ", chi_value,  ", df = ", df_value, ", P = ", p_value)}


create_data_structure <- function(branch_names, values, data_type = "branch") {
  n_branches <- length(branch_names)
  expected_length <- n_branches * 2 
  if (length(values) != expected_length) {
    warning(paste("values should be", expected_length, 
                  "but", length(values)))
    n_branches <- length(values) / 2
    if (n_branches != length(branch_names)) {
      branch_names <- paste0("Branch", 1:n_branches)}}
  
  data.frame(
    branch = factor(rep(branch_names, each = 2), levels = branch_names),
    group = factor(rep(c("early", "late"), n_branches), levels = c("early", "late")),number = values)}

# ==============================
# data
# ==============================
branch_names_3 <- c("Branch1", "Branch2", "Branch3")  
branch_names_4 <- c("Branch1", "Branch2", "Branch3", "Branch4")  
level_names <- c("Primary branch", "Secondary branch", "Tertiary branch") 

# Plant1 data
Plant1_branch_data <- c(3, 2, 2, 8, 3, 6)  
Plant1_branch_level_data <- c(0, 0, 7, 8, 1, 8) 
# Plant2 data
Plant2_branch_data <- c(5, 24, 4, 12, 1, 4, 1, 4)  
Plant2_branch_level_data <- c(1, 8, 7, 31, 3, 5)  
# Plant3 data
Plant3_branch_data <- c(3, 21, 2, 14, 1, 5)  
Plant3_branch_level_data <- c(2, 12, 2, 16, 2, 12)  

# Plant1 dataframe
Plant1_branch <- create_data_structure(branch_names_3, Plant1_branch_data)
Plant1_branch_level <- create_data_structure(level_names, Plant1_branch_level_data)
# Plant2 dataframe
Plant2_branch <- create_data_structure(branch_names_4, Plant2_branch_data)
Plant2_branch_level <- create_data_structure(level_names, Plant2_branch_level_data)
# Plant3 dataframe
Plant3_branch <- create_data_structure(branch_names_3, Plant3_branch_data)
Plant3_branch_level <- create_data_structure(level_names, Plant3_branch_level_data)

# ==============================
# Chi-square test
# ==============================

create_matrix_from_dataframe <- function(df) {
  n_branches <- length(unique(df$branch))
  early_values <- df$number[df$group == "early"]
  late_values <- df$number[df$group == "late"]
  matrix(c(early_values, late_values), nrow = n_branches, byrow = FALSE)}


chi_Plant1_branch <- calculate_chi_text(create_matrix_from_dataframe(Plant1_branch))
chi_Plant1_branch_level <- calculate_chi_text(create_matrix_from_dataframe(Plant1_branch_level))
chi_Plant2_branch <- calculate_chi_text(create_matrix_from_dataframe(Plant2_branch))
chi_Plant2_branch_level <- calculate_chi_text(create_matrix_from_dataframe(Plant2_branch_level))
chi_Plant3_branch <- calculate_chi_text(create_matrix_from_dataframe(Plant3_branch))
chi_Plant3_branch_level <- calculate_chi_text(create_matrix_from_dataframe(Plant3_branch_level))

# ==============================
# plot
# ==============================

plot_list <- list()

# Plant1 
plot_list$Plant1_branch <- create_plot(Plant1_branch, "plant1", chi_Plant1_branch, FALSE)
plot_list$Plant1_branch_level <- create_plot(Plant1_branch_level, "plant1", chi_Plant1_branch_level, FALSE)
# Plant2 
plot_list$Plant2_branch <- create_plot(Plant2_branch, "plant2", chi_Plant2_branch, FALSE)
plot_list$Plant2_branch_level <- create_plot(Plant2_branch_level, "plant2", chi_Plant2_branch_level, FALSE)
# Plant3 
plot_list$Plant3_branch <- create_plot(Plant3_branch, "plant3", chi_Plant3_branch, FALSE)
plot_list$Plant3_branch_level <- create_plot(Plant3_branch_level, "plant3", chi_Plant3_branch_level, FALSE)


for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
  cat("--- plot:", names(plot_list)[i], "---\n")
}

output_dir <- "../../"
dir.create(output_dir, showWarnings = FALSE)
