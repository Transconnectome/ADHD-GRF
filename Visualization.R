library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggpattern)
library(patchwork)


#### ATE plot ####
ate_df <- rbind(multi_env$ate, multi_env_gen$ate, multi_env_br$ate, multi_env_gen_br$ate,
                eur_env$ate, eur_env_gen$ate, eur_env_br$ate, eur_env_gen_br$ate)
ate_df$data <- as.factor(rep(c('Env','Env + Gen','Env + Br','Env + Gen + Br'), 2))
ate_df$sample <- as.factor(c(rep('Multiethnic sample', 4), rep('European sample', 4)))
#write.csv(ate_df,"ate.csv", row.names = FALSE)


data <- ate_df
# Fix column names (if needed)
colnames(data)[which(colnames(data) == "[95% ")] <- "lower_CI"
colnames(data)[which(colnames(data) == " CI]")] <- "upper_CI"

# Create the plot
ggplot(data, aes(x = sample, y = estimate, color = sample, alpha = data)) +
  geom_point(position = position_dodge(0.5), size = 4, show.legend = TRUE) +  # Point estimates with dodge for spacing
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                position = position_dodge(0.5), width = 0.2, size = 1.5, show.legend = FALSE) +  # Error bars
  scale_color_brewer(palette = "Dark2") +  # Use Dark2 color palette +
  scale_alpha_manual(values = c("Env" = 0.2, "Env + Br" = 0.4, "Env + Gen" = 0.6, "Env + Gen + Br" = 1)) +  # Adjust transparency for all four values
  ylim(0, 0.65) +
  labs(y = "ATE estimate", color = "Sample", alpha = "Model") +  # Customize labels
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        panel.grid.major.x = element_blank(),
        text = element_text(size = 16),
        legend.position = "none")





#### Test calibration plot ####
test_cal_df <- data.frame(rbind(multi_env$test_cal, multi_env_gen$test_cal, multi_env_br$test_cal, multi_env_gen_br$test_cal,
                     eur_env$test_cal, eur_env_gen$test_cal, eur_env_br$test_cal, eur_env_gen_br$test_cal))
test_cal_df$type <- rep(c('ate','hte'), 8)
test_cal_df$data <- as.factor(c(rep('Env',2),rep('Env + Gen',2),rep('Env + Br',2),rep('Env + Gen + Br',2)))
test_cal_df$sample <- as.factor(c(rep('Multiethnic sample', 8), rep('European sample', 8)))
test_cal_df_ate <- test_cal_df %>% filter(type == 'ate')
test_cal_df_hte <- test_cal_df %>% filter(type == 'hte')

#write.csv(test_cal_df_ate,"testcal_ate.csv", row.names = FALSE)
#write.csv(test_cal_df_hte,"testcal_hte.csv", row.names = FALSE)


data <- test_cal_df_hte
lower_CI <- data$Estimate - 1.96 * data$Std..Error
upper_CI <- data$Estimate + 1.96 * data$Std..Error
data <- data %>%
  mutate(significance = ifelse(Pr..t. < 0.05, "Significant", "Not Significant"))

# Create the HTE plot
data$sample <- factor(data$sample, levels = c("European sample", "Multiethnic sample"))  # Adjust sample order if needed
data$data <- factor(data$data, levels = c("Env", "Env + Br", "Env + Gen", "Env + Gen + Br"))  # Specify the order of levels

ggplot(data, aes(x = sample, y = Estimate, color = sample, alpha = data)) +
  geom_point(position = position_dodge(width = 0.5), size = 4, stroke = 2, show.legend = FALSE, shape = 16) +  # Points with dodge
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, position = position_dodge(width = 0.5), linewidth = 1.5, show.legend = FALSE) +  # Error bars with dodge
  #scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +  # Shape for significance
  scale_alpha_manual(values = c("Env" = 0.2, "Env + Br" = 0.4, "Env + Gen" = 0.6, "Env + Gen + Br" = 1)) +  # Alpha transparency
  scale_color_brewer(palette = "Dark2") +  # Use Dark2 color palette +
  #scale_color_manual(values = c("European" = "#1B9E77", "Multiethnic" = "#D95F02")) +  # Colors for samples
  labs(x = "", y = "ITE coefficient") +  # Labels
  scale_y_continuous(breaks = seq(-7.5, 5, by = 2.5)) +  # Y-axis breaks
  theme_minimal() +  # Minimal theme
  theme(panel.grid.major.x = element_blank(),
        text = element_text(size = 16),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")  # Remove the legend


  
# Create the ATE plot
ggplot(data, aes(x = sample, y = Estimate, shape = significance, color = data)) +
  geom_point(position = position_dodge(width = 0.5), size = 4, stroke = 2) +  # Add points with dodge for separation
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, position = position_dodge(width = 0.5), , linewidth = 1.5) +  # Add error bars with dodge
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +  # Filled for significant, empty for non-significant
  labs(x = "", y = "ATE coefficient") +  # Set axis labels
  scale_color_brewer(palette = "Dark2") +  # Use the Dark2 palette for color-coding
  ylim(0, 2) +
  theme_minimal() +  # Use a clean theme
  theme(legend.title = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 16),
        panel.grid.minor.y = element_blank()) + 
  guides(color = guide_legend(title = "Model"))  # Add a legend for the data categories



#### Variable importance plot ####
var_imp_multi <- multi_env_gen$var_imp
var_imp_eur <- eur_env_gen$var_imp
var_imp_df <- var_imp_multi
names(var_imp_df) <- 'multi'
var_imp_df$eur <- 0
var_imp_df$group <- NA 
var_imp_df$group[which(rownames(var_imp_df) %in% c('interview_age','sex','demo_comb_income_v2',
                                                 'ethnic_g_AFR','ethnic_g_AMR','ethnic_g_EAS','ethnic_g_EUR'))] <- 'demographic'
var_imp_df$group[which(rownames(var_imp_df) %in% names(PRS)[2:31])] <- 'genetic'
var_imp_df$group[which(is.na(var_imp_df$group))] <- 'environmental'

for(i in 1:nrow(var_imp_eur)){
  ind <- which(rownames(var_imp_multi) == rownames(var_imp_eur)[i])
  var_imp_df$eur[ind] <- var_imp_eur$importance[i] 
  
  if(rownames(var_imp_eur)[i] %in% c('BMIeur4','DEPeur4','ALCDEP_EURauto','PTSDeur4','SCZ_EURauto')){
    ind <- which(substr(rownames(var_imp_multi), 1, 3) == substr(rownames(var_imp_eur)[i], 1, 3))
    var_imp_df$eur[ind] <- var_imp_eur$importance[i] 
  }
}


#write.csv(var_imp_df,"var_imp.csv", row.names = T)

data <- var_imp_df
data$Variable <- rownames(data)

# Sort the data by the sum of Multiethnic and European values (optional, for better visualization)
data <- data %>%
  arrange(desc(multi + eur))
data <- data[1:20,]

# Assuming your dataframe is called data
data_long <- data %>%
  pivot_longer(cols = c(multi, eur),  # Specify the columns to merge
               names_to = "sample",   # New column to indicate source (multi/eur)
               values_to = "imp")     # New column for the merged values


library(ggpattern)

ggplot(data) +
  # Plot Multiethnic sample (negative values for left side) with pattern for group
  geom_bar_pattern(aes(x = reorder(Variable, multi + eur, decreasing = FALSE),
                       y = -eur, pattern = group, alpha = 0.9), 
                   fill = "#1B9E77", stat = "identity", 
                   pattern_fill = "#1B9E77", 
                   pattern_density = 0.3, 
                   pattern_spacing = 0.05,
                   pattern_angle = 30,
                   pattern_colour = NA) +
  
  # Plot European sample (positive values for right side) with pattern for group
  geom_bar_pattern(aes(x = Variable, y = multi, pattern = group, alpha = 0.9), 
                   fill = "#D95F02", stat = "identity", 
                   pattern_fill = "#D95F02", 
                   pattern_density = 0.3, 
                   pattern_spacing = 0.05,
                   pattern_angle = 30,
                   pattern_colour = NA) +
  
  # Add horizontal threshold lines
  geom_hline(yintercept = -eur_env_gen$threshold, color = "black", linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = multi_env_gen$threshold, color = "black", linetype = "dashed", size = 0.8) +
  
  # Flip coordinates for better readability
  coord_flip() +
  
  # Add a vertical line at x = 0 to separate the two samples
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.8) +
  
  # Assign specific patterns for each group using scale_pattern_manual()
  scale_pattern_manual(values = c("environmental" = "stripe", 
                                  "genetic" = "none", 
                                  "Group3" = "polka dot", 
                                  "Group4" = "circle")) +
  
  # Labels and theme adjustments
  labs(x = NULL, y = "Variable Importance", fill = "Variable Types", title = "") +
  
  # Adjust theme and remove unnecessary elements
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),  # Remove y-axis ticks
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")  # Adjust legend position



#### Longitudinal outcomes ATE plot ####
ate_df <- rbind(multi_env_gen_1y$ate, multi_env_gen_2y$ate, multi_env_gen_3y$ate)
ate_df$data <- as.factor(c('1 year','2 year','3 year'))


data <- ate_df
# Fix column names (if needed)
colnames(data)[which(colnames(data) == "[95% ")] <- "lower_CI"
colnames(data)[which(colnames(data) == " CI]")] <- "upper_CI"
color <- c("1 year" = "#66a61e", "2 year" = "#e6ab02", "3 year" = "#a6761d")

# Create the plot
ggplot(data, aes(x = data, y = estimate, alpha = data)) +
  geom_point(size = 4, color = '#7570B3') +  # Set a constant purple color for points
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, size = 1.5, color = '#7570B3') +  # Set a constant purple color for error bars
  scale_x_discrete(expand = c(0.5, 0.05)) +
  scale_alpha_manual(values = c("1 year" = 0.3, "2 year" = 0.6, "3 year" = 1.0)) +
  ylim(0, 0.65) +
  labs(y = "ATE estimate", color = "Model") +  # Customize labels
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(size = 16))

#write.csv(ate_df,"ate_long.csv", row.names = FALSE)



#### Test calibration plot ####
test_cal_df <- data.frame(rbind(multi_env_gen_1y$test_cal, multi_env_gen_2y$test_cal, multi_env_gen_3y$test_cal))
test_cal_df$type <- rep(c('ate','hte'), 3)
test_cal_df$data <- as.factor(c(rep('1 year',2),rep('2 year',2),rep('3 year',2)))
test_cal_df_ate <- test_cal_df %>% filter(type == 'ate')
test_cal_df_hte <- test_cal_df %>% filter(type == 'hte')

#write.csv(test_cal_df_ate,"testcal_ate_long.csv", row.names = FALSE)
#write.csv(test_cal_df_hte,"testcal_hte_long.csv", row.names = FALSE)


data <- test_cal_df_hte
lower_CI <- data$Estimate - 1.96 * data$Std..Error
upper_CI <- data$Estimate + 1.96 * data$Std..Error
data <- data %>%
  mutate(significance = ifelse(Pr..t. < 0.05, "Significant", "Not Significant"))

# Create the HTE plot
ggplot(data, aes(x = data, y = Estimate, alpha = data)) +
  geom_point(size = 4, stroke = 2, color = '#7570B3', shape=16) +  # Add points with dodge for separation
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, linewidth = 1.5, color = '#7570B3') +  # Add error bars with dodge
  scale_x_discrete(expand = c(0.5, 0.05)) +
  scale_alpha_manual(values = c("1 year" = 0.3, "2 year" = 0.6, "3 year" = 1.0)) +
  #scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +  # Filled for significant, empty for non-significant
  labs(x = "", y = "HTE coefficient") +  # Set axis labels
  #scale_color_manual(values = color) +
  scale_y_continuous(breaks = seq(-10, 5, by = 2.5)) + 
  ylim(c(-10,5)) +
  theme_minimal() +  # Use a clean theme
  theme(legend.title = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 16),
        panel.grid.minor.y = element_blank()) + 
  #geom_hline(yintercept = 0, color = "black", linetype = "solid")
guides(color = guide_legend(title = "Model"))  # Add a legend for the data categories




#### partial dependence ####
univ.sim <- rbind(univ.sim_eur, univ.sim_multi, univ.sim_multi_1y)
univ.sim$sample <- factor(c(rep('European', nrow(univ.sim_eur)), 
                          rep('Multiethnic',nrow(univ.sim_multi)), 
                          rep('Multiethnic 1-year outcome', nrow(univ.sim_multi_1y))))

ind <- which(names(univ.sim) %in% c('trend'))
univ.sim <- univ.sim[,-ind] 

#[1] "sleep_prob"   "SMOKERauto"   "RISKTOLauto"  "BMImulti"     "SCZ_METAauto" "WORRYauto"   
#[7] "EDauto"      

select_feature <- 'EDauto'
univ_sim_data <- univ.sim %>% filter(feature == select_feature)
univ_sim_data_gather <- data.frame(matrix(ncol = 4))
univ_sim_data_gather <- univ_sim_data_gather[-1,]
colnames(univ_sim_data_gather) <- c('feature','sample','quantile','estimate')

for(ind in 1:nrow(univ_sim_data)){
  add <- univ_sim_data[ind,] %>%
    gather(quantile, estimate, tau.hat.Q0:tau.hat.Q4)
  add$quantile <- c('Q0','Q1','Q2','Q3','Q4')
  univ_sim_data_gather <- rbind(univ_sim_data_gather, add)
}

color <- c('European' = '#1b9e77', 'Multiethnic' = '#d95f02', 'Multiethnic 1-year outcome' = '#7570b3')

ggplot(univ_sim_data_gather, aes(x=quantile, y=estimate, group = sample, color = sample)) +
  geom_point(size=5, show.legend = FALSE) + 
  geom_line(linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = color) +
  ylim(c(0.3, 0.8)) +
  #scale_y_continuous(limits=c(0.3, 0.5)) +
  #ggtitle(univ_sim_data$feature[1]) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_blank())



#### Longitudinal variable importance plot ####
var_imp_df <- multi_env_gen_1y$var_imp
names(var_imp_df) <- 'multi'
var_imp_df$group <- NA 
var_imp_df$group[which(rownames(var_imp_df) %in% c('interview_age','sex','demo_comb_income_v2',
                                                   'ethnic_g_AFR','ethnic_g_AMR','ethnic_g_EAS','ethnic_g_EUR'))] <- 'demographic'
var_imp_df$group[which(rownames(var_imp_df) %in% names(PRS)[2:31])] <- 'genetic'
var_imp_df$group[which(is.na(var_imp_df$group))] <- 'environmental'

#write.csv(var_imp_df,"var_imp_long.csv", row.names = T)

data <- var_imp_df
data$Variable <- rownames(data)
data <- data[1:20,]
library(ggpattern)

ggplot(data) +
  # Plot Multiethnic sample with patterns for the "group" variable
  geom_bar_pattern(aes(x = reorder(Variable, multi, decreasing = FALSE), 
                       y = multi, pattern = group), 
                   fill = "#7570B3", alpha = 0.7, stat = "identity",
                   pattern_density = 0.3, pattern_spacing = 0.05, 
                   pattern_fill = "#7570B3", pattern_color = NA,
                   pattern_angle = 30) +

  # Add threshold line for multiethnic sample
  geom_hline(yintercept = multi_env_gen_1y$threshold, color = "black", linetype = "dashed", size = 0.8) +
  
  # Flip coordinates for better readability
  coord_flip() +
  
  # Remove the color coding for group, apply patterns instead
  scale_fill_manual(values = c("#7570B3")) +  # All bars in purple
  
  scale_pattern_manual(values = c("environmental" = "stripe", 
                                  "genetic" = "none", 
                                  "demographic" = "circle")) +
  
  # Add vertical line at x = 0
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.8) +
  
  # Labels and theme adjustments
  labs(x = NULL, y = "Variable Importance", fill = "Variable Types", title = "") +
  
  # Adjust theme and remove unnecessary elements
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),  # Remove y-axis ticks
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right")  # Adjust legend position


############################## @ Figure 1 #######################################
library(ggplot2)
#detach("package:ggstatsplot", unload = TRUE)
library(patchwork)
library(grid)
library(dplyr)

#========== Plot 1 ==========
ate_df <- rbind(multi_env$ate, multi_env_gen$ate, multi_env_br$ate, multi_env_gen_br$ate,
                eur_env$ate, eur_env_gen$ate, eur_env_br$ate, eur_env_gen_br$ate)
ate_df$data <- as.factor(rep(c('Env','Env + Gen','Env + Br','Env + Gen + Br'), 2))
ate_df$sample <- as.factor(c(rep('Multiethnic', 4), rep('European', 4)))
ate_df <- ate_df %>%
  mutate(stars = case_when(
    `P-value` < 0.001 ~ "***",
    `P-value` < 0.01 ~ "**",
    `P-value` < 0.05 ~ "*",
    TRUE ~ ""
  ))

data <- ate_df
colnames(data)[which(colnames(data) == "[95% ")] <- "lower_CI"
colnames(data)[which(colnames(data) == " CI]")] <- "upper_CI"

p1 <- ggplot(data, aes(x = sample, y = estimate, color = sample, alpha = data)) +
  geom_point(position = position_dodge(0.8), size = 2.5, shape = 16, show.legend = TRUE) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
                position = position_dodge(0.8), width = 0.4, linewidth = 0.6, show.legend = FALSE) +
  geom_text(aes(label = stars, y = upper_CI + 0.01),  # small space above error bar
            position = position_dodge(0.8),
            size = 3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_alpha_manual(values = c("Env" = 0.5, "Env + Br" = 0.6, "Env + Gen" = 0.8, "Env + Gen + Br" = 1)) +
  ylim(0, 0.65) +
  labs(title = "Average effect of parental ADHD \non baseline child ADHD",
       y = "Average effect of parental ADHD (ATE)",
       color = "Sample", alpha = "Model") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(family = "Arial", size = 7),
    plot.title = element_text(hjust = 0.5, size = 7),
    legend.position = "none",
  )

# ========== Plot 2 ==========
test_cal_df <- data.frame(rbind(multi_env$test_cal, multi_env_gen$test_cal, multi_env_br$test_cal, multi_env_gen_br$test_cal,
                                eur_env$test_cal, eur_env_gen$test_cal, eur_env_br$test_cal, eur_env_gen_br$test_cal))
test_cal_df$type <- rep(c('ate','hte'), 8)
test_cal_df$data <- as.factor(c(rep('Env',2),rep('Env + Gen',2),rep('Env + Br',2),rep('Env + Gen + Br',2)))
test_cal_df$sample <- as.factor(c(rep('Multiethnic', 8), rep('European', 8)))
test_cal_df_ate <- test_cal_df %>% filter(type == 'ate')
test_cal_df_hte <- test_cal_df %>% filter(type == 'hte')
test_cal_df_hte <- test_cal_df_hte %>%
  mutate(stars = case_when(
    `Pr..t.` < 0.001 ~ "***",
    `Pr..t.` < 0.01 ~ "**",
    `Pr..t.` < 0.05 ~ "*",
    TRUE ~ ""
  ))


data <- test_cal_df_hte
data <- data %>%
  mutate(
    lower_CI = Estimate - 1.96 * Std..Error,
    upper_CI = Estimate + 1.96 * Std..Error,
    significance = ifelse(Pr..t. < 0.05, "Significant", "Not Significant")
  )

p2 <- ggplot(data, aes(x = sample, y = Estimate, color = sample, alpha = data)) +
  geom_point(position = position_dodge(0.8), size = 2.5, shape = 16, show.legend = FALSE) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
                width = 0.4, position = position_dodge(0.8), linewidth = 0.6, show.legend = FALSE) +
  geom_text(aes(label = stars, y = upper_CI + 0.15),  # small space above error bar
            position = position_dodge(0.8),
            size = 3, show.legend = FALSE) +  
  scale_alpha_manual(values = c("Env" = 0.5, "Env + Br" = 0.6, "Env + Gen" = 0.8, "Env + Gen + Br" = 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Estimated significance of individual-level effect \nof parental ADHD on baseline child ADHD",
       y = "Individual effect coefficient (CATE coefficient)") +
  #scale_y_continuous(breaks = seq(-15.0, 5, by = 2.5)) +
  ylim(-15, 5) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, size = 7),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

#==================== Plot 3-1 ====================
ate_df <- rbind(multi_env_gen_1y$ate, multi_env_gen_2y$ate, multi_env_gen_3y$ate)
ate_df$data <- as.factor(c('1 year','2 year','3 year'))
colnames(ate_df)[which(colnames(ate_df) == "[95% ")] <- "lower_CI"
colnames(ate_df)[which(colnames(ate_df) == " CI]")] <- "upper_CI"
ate_df <- ate_df %>%
  mutate(stars = case_when(
    `P-value` < 0.001 ~ "***",
    `P-value` < 0.01 ~ "**",
    `P-value` < 0.05 ~ "*",
    TRUE ~ ""
  ))

p3_1 <- ggplot(ate_df, aes(x = data, y = estimate, alpha = data)) +
  geom_point(size = 2.5, stroke = 0.5, color = '#7570B3', shape = 16) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.4, linewidth = 0.6, color = '#7570B3') +
  geom_text(aes(label = stars, y = upper_CI + 0.01),  # small space above error bar
            position = position_dodge(0.8),
            size = 3, show.legend = FALSE, color = '#7570B3') +  
  scale_x_discrete(expand = c(0.5, 0.05)) +
  scale_alpha_manual(values = c("1 year" = 0.6, "2 year" = 0.8, "3 year" = 1.0)) +
  ylim(0, 0.65) +
  labs(y = "Average effect of parental ADHD (ATE)") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(size = 7),
    legend.position = "none",
    axis.text.x = element_blank()
  )

#==================== Plot 3-2 ====================
test_cal_df <- data.frame(rbind(multi_env_gen_1y$test_cal, multi_env_gen_2y$test_cal, multi_env_gen_3y$test_cal))
test_cal_df$type <- rep(c('ate','hte'), 3)
test_cal_df$data <- as.factor(c(rep('1 year',2),rep('2 year',2),rep('3 year',2)))
test_cal_df_hte <- test_cal_df %>% filter(type == 'hte') %>%
  mutate(
    lower_CI = Estimate - 1.96 * Std..Error,
    upper_CI = Estimate + 1.96 * Std..Error,
    significance = ifelse(Pr..t. < 0.05, "Significant", "Not Significant")
  )
test_cal_df_hte <- test_cal_df_hte %>%
  mutate(stars = case_when(
    `Pr..t.` < 0.001 ~ "***",
    `Pr..t.` < 0.01 ~ "**",
    `Pr..t.` < 0.05 ~ "*",
    TRUE ~ ""
  ))

p3_2 <- ggplot(test_cal_df_hte, aes(x = data, y = Estimate, alpha = data)) +
  geom_point(size = 2.5, stroke = 0.5, color = '#7570B3', shape = 16) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.4, linewidth = 0.6, color = '#7570B3') +
  geom_text(aes(label = stars, y = upper_CI + 0.15),  # small space above error bar
            position = position_dodge(0.8),
            size = 3, show.legend = FALSE, color = '#7570B3') +  
  scale_x_discrete(expand = c(0.5, 0.05)) +
  scale_alpha_manual(values = c("1 year" = 0.6, "2 year" = 0.8, "3 year" = 1.0)) +
  ylim(-15, 5) +
  #scale_y_continuous(breaks = seq(-10, 5, 2.5)) +
  labs(y = "Individual effect coefficient (CATE coefficient)") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(size = 7),
    legend.position = "none",
    axis.text.x = element_blank()
  )

p3 <- (p3_1 + p3_2) +
  plot_layout(nrow = 1) +  # Ensure layout is defined for subplots
  plot_annotation(
    title = "Average and individual-level effect of parental ADHD\non longitudinal child ADHD in the multiethnic sample",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 7, family = "Arial"),
      plot.margin = margin(5, 5, 5, 5)  # optional, keeps spacing clean
    )
  )


#==================== Combine and Save ====================
final_plot <- p1 + p2 + p3 + plot_layout(nrow = 1, widths = c(1, 1, 1.2))

ggsave("ADHD_combined_figur31e.pdf", final_plot,
       device = cairo_pdf, width = 7.08, height = 3.5, units = "in", dpi = 600)


