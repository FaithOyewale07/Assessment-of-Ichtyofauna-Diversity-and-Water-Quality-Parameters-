# ==========================================
# ALPHA DIVERSITY INDICES WITH POOLED
# Including All Stations Combined
# ==========================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(vegan)
library(tidyverse)
library(grid)

# Calculate diversity for POOLED data (all stations combined)
pooled_data <- colSums(plankton_data_new)  # Sum across all stations
pooled_matrix <- matrix(pooled_data, nrow = 1)  # Convert to matrix format
colnames(pooled_matrix) <- names(pooled_data)

# Calculate pooled diversity indices
diversity_pooled <- data.frame(
  stations = "Pooled",
  Shannon = round(diversity(pooled_matrix, index = "shannon"), 2),
  Simpson = round(diversity(pooled_matrix, index = "simpson"), 2),
  InvSimpson = round(diversity(pooled_matrix, index = "invsimpson"), 2),
  Richness = specnumber(pooled_matrix),
  Pielou_Evenness = round(diversity(pooled_matrix, index = "shannon") / 
    log(specnumber(pooled_matrix)), 2),
  Margalef = round((specnumber(pooled_matrix) - 1) / log(sum(pooled_matrix) + 1),2),
  Total_Abundance = sum(pooled_matrix)
)

print("Pooled Diversity Indices:")
print(diversity_pooled)

# Combine individual and pooled
diversity_with_pooled <- bind_rows(diversity_indices, diversity_pooled) |> 
  mutate(stations = factor(stations, levels = c("Ibeshe", "Ilashe", "Igboalejo", "Igboeseyero", 
  "Igboologun", "Pooled")))

print("\nDiversity Indices with Pooled:")
print(diversity_with_pooled)

write.csv(diversity_with_pooled, 'diversity_complete.csv')

# ==========================================
# STEP 2: RESHAPE DATA FOR PLOTTING
# ==========================================

diversity_plot_data <- diversity_with_pooled |> 
  pivot_longer(cols = c(Shannon, Simpson, InvSimpson, Richness, 
                        Pielou_Evenness, Margalef, Total_Abundance),
               names_to = "Index",
               values_to = "Value")

# Define colors (5 stations + 1 pooled)
station_colors <- c("Ibeshe" = "#4DBBD5", 
                    "Ilashe" = "#E64B35",
                    "Igboalejo" = "#00A087",
                    "Igboeseyero" = "#F39B7F",
                    "Igboologun" = "#8491B4",
                    "Pooled" = "#7E6148")  # Brown color for pooled

# ==========================================
# STEP 3: STATISTICAL TESTS (EXCLUDING POOLED)
# ==========================================

# Calculate stats only for individual stations (not pooled)
diversity_no_pooled <- diversity_plot_data |> 
  filter(stations != "Pooled")

stat_results <- diversity_no_pooled |> 
  group_by(Index) |> 
  summarise(
    p_value = kruskal.test(Value ~ stations)$p.value,
    p_label = case_when(
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01 ~ paste0("p = ", sprintf("%.3f", p_value)),
      p_value < 0.05 ~ paste0("p = ", sprintf("%.2f", p_value)),
      TRUE ~ paste0("p = ", sprintf("%.2f", p_value))
    )
  )

print("\nStatistical Results (Kruskal-Wallis - excluding Pooled):")
print(stat_results)

write.csv(stat_results, 'statistic_result_no_pooled.csv')

# ==========================================
# STEP 4:  DIVERSITY WITH POOLED
# ==========================================

cat("\n=== CREATING PANEL A WITH POOLED ===\n")

# Create individual plots
plot_list <- list()
index_order <- c("Shannon", "Simpson", "Richness", "Pielou_Evenness",
                 "Margalef", "InvSimpson", "Total_Abundance")

for(idx in index_order) {
  
  data_subset <- diversity_plot_data %>% filter(Index == idx)
  p_val <- stat_results %>% filter(Index == idx) %>% pull(p_label)
  
  # Highlight pooled with different alpha
  p <- ggplot(data_subset, aes(x = stations, y = Value, fill = stations)) +
    geom_col(aes(alpha = stations == "Pooled"), 
             width = 0.6, color = "black", size = 0.7) +
    scale_alpha_manual(values = c("TRUE" = 0.95, "FALSE" = 0.75), guide = "none") +
    scale_fill_manual(values = station_colors) +
    theme_bw() +
    labs(title = idx, y = "Value", x = NULL) +
    annotate("text", x = 3, y = max(data_subset$Value, na.rm = TRUE) * 0.92, 
             label = p_val, size = 3.2, fontface = "italic") +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  plot_list[[idx]] <- p
}

# Arrange in 2x4 grid
p_diversity_panel <- grid.arrange(
  grobs = plot_list,
  ncol = 4,
  nrow = 2,
  top = textGrob("A  Alpha Diversity Indices (Individual Stations + Pooled)", 
                 gp = gpar(fontface = "bold", fontsize = 14),
                 hjust = 0, x = 0.02)
)

ggsave("PanelA_Diversity_with_Pooled.png", 
       plot = p_diversity_panel, 
       width = 16, height = 8, dpi = 300)

# ==========================================
# STEP 5: ABUNDANT TAXA WITH POOLED
# ==========================================

cat("\n\n=== ABUNDANT TAXA WITH POOLED ===\n")

# Use same classification as before
total_abundance_per_taxon <- colSums(plankton_data_new)
relative_abundance_per_taxon <- round((total_abundance_per_taxon / sum(total_abundance_per_taxon)) * 100, 2)

taxa_classification <- data.frame(
  Taxa = names(total_abundance_per_taxon),
  Total_Abundance = total_abundance_per_taxon,
  Relative_Abundance_Percent = relative_abundance_per_taxon
) |> 
  mutate(
    Classification = case_when(
      Relative_Abundance_Percent >= 5.0 ~ "Abundant",
      Relative_Abundance_Percent <= 0.5 ~ "Rare",
      TRUE ~ "Intermediate"
    )
  ) |> 
  arrange(desc(Total_Abundance))

# Summary
classification_summary <- taxa_classification |> 
  group_by(Classification) |> 
  summarise(
    Number_of_Taxa = n(),
    Total_Abundance = round(sum(Total_Abundance),2),
    Mean_Rel_Abundance = round(mean(Relative_Abundance_Percent), 2),
    Percent_of_Total = round(sum(Total_Abundance) / sum(taxa_classification$Total_Abundance) * 100, 2),
  )
print(classification_summary)

# Extract abundant and rare taxa
abundant_data <- plankton_data_new |> 
  select(all_of(abundant_taxa))

abundant_taxa <- taxa_classification |>  
  filter(Classification == "Abundant") |>  
  pull(Taxa)

intermediate_taxa <- taxa_classification |> 
  filter(Classification == 'Intermediate') |> 
  pull(Taxa)

rare_taxa <- taxa_classification |> 
  filter(Classification == "Rare") |> 
  pull(Taxa)

# Individual stations
abundant_data_individual <- plankton_data_new[, abundant_taxa, drop = FALSE]

abundant_diversity_individual <- data.frame(
  stations = metadata$stations,
  Shannon = diversity(abundant_data_individual, index = "shannon"),
  Simpson = diversity(abundant_data_individual, index = "simpson"),
  InvSimpson = diversity(abundant_data_individual, index = "invsimpson"),
  Richness = specnumber(abundant_data_individual),
  Evenness = diversity(abundant_data_individual, index = "shannon") / 
    log(specnumber(abundant_data_individual)),
  Total_Abundance = rowSums(abundant_data_individual)
)

# Pooled abundant taxa
abundant_pooled_data <- colSums(abundant_data_individual)
abundant_pooled_matrix <- matrix(abundant_pooled_data, nrow = 1)
colnames(abundant_pooled_matrix) <- names(abundant_pooled_data)

abundant_diversity_pooled <- data.frame(
  stations = "Pooled",
  Shannon = diversity(abundant_pooled_matrix, index = "shannon"),
  Simpson = diversity(abundant_pooled_matrix, index = "simpson"),
  InvSimpson = diversity(abundant_pooled_matrix, index = "invsimpson"),
  Richness = specnumber(abundant_pooled_matrix),
  Evenness = diversity(abundant_pooled_matrix, index = "shannon") / 
    log(specnumber(abundant_pooled_matrix)),
  Total_Abundance = sum(abundant_pooled_matrix)
)

# Combine
abundant_diversity_with_pooled <- bind_rows(
  abundant_diversity_individual, 
  abundant_diversity_pooled
) |> 
  mutate(stations = factor(stations, levels = c("Ibeshe", "Ilashe", "Igboalejo", "Igboeseyero", "Igboologun", "Pooled")))

# Replace NaN with 0
abundant_diversity_with_pooled[is.na(abundant_diversity_with_pooled)] <- 0

print("Abundant Taxa Diversity with Pooled:")
print(abundant_diversity_with_pooled)

# Reshape
abundant_long <- abundant_diversity_with_pooled |> 
  pivot_longer(cols = c(Shannon, Simpson, InvSimpson, Richness, Evenness, Total_Abundance),
               names_to = "Index",
               values_to = "Value")

# Statistical tests (excluding pooled)
abundant_no_pooled <- abundant_long |> 
  filter(stations != "Pooled")

abundant_stats <- abundant_no_pooled |> 
  group_by(Index) |> 
  summarise(
    p_value = kruskal.test(Value ~ stations)$p.value,
    p_label = case_when(
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01 ~ paste0("p = ", sprintf("%.3f", p_value)),
      TRUE ~ paste0("p = ", sprintf("%.2f", p_value))
    )
  )

# Individual abundant taxa abundance by station
abundant_taxa_long <- plankton_data_new |> 
  select(all_of(abundant_taxa)) |> 
  mutate(stations = metadata$stations) |> 
  pivot_longer(cols = all_of(abundant_taxa),
               names_to = "Taxa",
               values_to = "Abundance")

# Plot individual abundant taxa
p_abundant_individual <- ggplot(abundant_taxa_long, 
                                aes(x = stations, y = Abundance, fill = stations)) +
  geom_col(alpha = 0.7, width = 0.5, color = "black", size = 0.4) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = station_colors) +
  theme_bw() +
  labs(title = "Individual Abundant Taxa Distribution (>5% Relative Abundance)",
       subtitle = paste("Abundance across Stations based on", length(abundant_taxa), "abundant taxa"),
       x = "Stations",
       y = "Abundance (cells/L)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 9, face = "italic"),
        strip.background = element_rect(fill = "lightblue"),
        legend.position = "bottom")

print(p_abundant_individual)

ggsave("adundant_Individual_taxa_.png", 
       plot = p_abundant_individual, width = 14, height = 9, dpi = 300)

#################################
## Intermediate taxa
#################################

# Individual stations
intermediate_data_individual <- plankton_data_new |> 
  select(all_of(intermediate_taxa))

intermediate_diversity_individual <- data.frame(
  stations = metadata$stations,
  Shannon = diversity(intermediate_data_individual, index = "shannon"),
  Simpson = diversity(intermediate_data_individual, index = "simpson"),
  InvSimpson = diversity(intermediate_data_individual, index = "invsimpson"),
  Richness = specnumber(intermediate_data_individual),
  Evenness = diversity(intermediate_data_individual, index = "shannon") / 
    log(specnumber(intermediate_data_individual)),
  Total_Abundance = rowSums(intermediate_data_individual)
)

# Pooled abundant taxa
intermediate_pooled_data <- colSums(intermediate_data_individual)
intermediate_pooled_matrix <- matrix(intermediate_pooled_data, nrow = 1)
colnames(intermediate_pooled_matrix) <- names(intermediate_pooled_data)

intermediate_diversity_pooled <- data.frame(
  stations = "Pooled",
  Shannon = diversity(intermediate_pooled_matrix, index = "shannon"),
  Simpson = diversity(intermediate_pooled_matrix, index = "simpson"),
  InvSimpson = diversity(intermediate_pooled_matrix, index = "invsimpson"),
  Richness = specnumber(intermediate_pooled_matrix),
  Evenness = diversity(intermediate_pooled_matrix, index = "shannon") / 
    log(specnumber(intermediate_pooled_matrix)),
  Total_Abundance = sum(intermediate_pooled_matrix)
)

# Combine
intermediate_diversity_with_pooled <- bind_rows(
  intermediate_diversity_individual, 
  intermediate_diversity_pooled
) |> 
  mutate(stations = factor(stations, levels = c("Ibeshe", "Ilashe", "Igboalejo", "Igboeseyero", "Igboologun", "Pooled")))

print("Intermediate Taxa Diversity with Pooled:")
print(Intermediate_diversity_with_pooled)

# Reshape
intermediate_long <- intermediate_diversity_with_pooled |> 
  pivot_longer(cols = c(Shannon, Simpson, InvSimpson, Richness, Evenness, Total_Abundance),
               names_to = "Index",
               values_to = "Value")

# Statistical tests (excluding pooled)
intermediate_no_pooled <- intermediate_long |> 
  filter(stations != "Pooled")

intermediate_stats <- intermediate_no_pooled |> 
  group_by(Index) |> 
  summarise(
    p_value = kruskal.test(Value ~ stations)$p.value,
    p_label = case_when(
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01 ~ paste0("p = ", sprintf("%.3f", p_value)),
      TRUE ~ paste0("p = ", sprintf("%.2f", p_value))
    )
  )

# Individual Intermediate taxa abundance by station
intermediate_taxa_long <- plankton_data_new |> 
  select(all_of(intermediate_taxa)) |> 
  mutate(stations = metadata$stations) |> 
  pivot_longer(cols = all_of(intermediate_taxa),
               names_to = "Taxa",
               values_to = "Abundance")

# Plot individual abundant taxa
p_intermediate_individual <- ggplot(intermediate_taxa_long$Taxa[1], 
                                aes(x = stations, y = Abundance, fill = stations)) +
  geom_col(alpha = 0.7, width = 0.5, color = "black", size = 0.4) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = station_colors) +
  theme_bw() +
  labs(title = "Individual intermediate Taxa Distribution (>0.5:<5% Relative Abundance)",
       subtitle = paste("Abundance across Stations based on", length(intermediate_taxa), "intermediate taxa"),
       x = "Stations",
       y = "Abundance (cells/L)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 9, face = "italic"),
        strip.background = element_rect(fill = "lightblue"),
        legend.position = "bottom")

print(p_intermediate_individual)

ggsave("intermediate_Individual_taxa.png", 
       plot = p_intermediate_individual, width = 14, height = 9, dpi = 300)

write.csv(intermediate_data_individual, 'intermediate_data_individual.csv')
write.csv(abundant_diversity_individual, 'abundant_diversity_individual.csv')
write.csv(intermediate_diversity_with_pooled, 'intermediate_diversity_with_pooled.csv')
write.csv(taxa_classification, 'taxa_classification.csv')
write.csv(intermediate_stats, 'intermediate_stats.csv')

