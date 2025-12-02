## Data Visualization of Alpha-Diversity

div_long <- diversity_indices |> 
  select(stations, Shannon, Simpson, Richness, Pielou_Evenness) |> 
  pivot_longer(
    cols = c(Shannon, Simpson, Richness, Pielou_Evenness),
    names_to = 'Variables',
    values_to = 'Value'
  )

## Plot Diversity Indices
plankton_diversity <- div_long |> 
  ggplot(aes(stations, Value, fill = Variables))+
  geom_bar(stat = "identity", position = 'dodge',
           width = 0.4, size = 0.25, color = "black") +
  facet_wrap(~Variables, scales = "free_y", ncol = 2) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")+
  labs(title = "Phytoplankton Diversity Indices by Station",
       y = "Index Value", x = "Station") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

print(plankton_diversity)

# Create boxplot with facets
p_boxplot_facet <- ggplot(div_long, aes(x = stations, y = Value, fill = stations)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21,  outlier.size = 2) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "white", color = "black") +
  facet_wrap(~Variables, ncol = 2) +
  theme_bw() +
  labs(title = "Phytoplankton Diversity Indices by Station",
       subtitle = "Boxplots with Individual Data Points",
       x = "Station",
       y = "Index Value") +
  scale_fill_brewer(palette = "Set3") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "lightgray"),
        legend.position = "none")

print(p_boxplot_facet)


p_shannon <- ggplot(diversity_indices, aes(x = stations, y = Shannon, fill = stations)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 4, alpha = 0.9, shape = 21, color = "black") +
  theme_bw() +
  labs(title = "Shannon Diversity Index (H')",
       x = "Station",
       y = "Shannon Index Value") +
  scale_fill_brewer(palette = "Set2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  ylim(0, max(diversity_indices$Shannon) * 1.1)

print(p_shannon)
