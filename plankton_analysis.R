## Load the required Libraries for Analysis
## Alpha Diversity Indices
library(tidyverse)
library(dplyr)
library(vegan)
library(reshape2)
library(readxl)
library(janitor)
library(car)
library(FSA)
library(psych)
library(betapart)
library(indicspecies)
library(devtools)
library(pairwiseAdonis)
library(ggrepel)
library(knitr)

### Import the data set 
plankton_data <- read_excel('PLANKTON BOSUN.xlsx', sheet = 2)

## Convert column names to snake case
plankton_data <- plankton_data |> 
  janitor::clean_names()
  

# Clean the data - remove division/class rows, keep only species
# Method 1: Keep only rows with species names 
## (rows that don't start with "Division:" or "Class:")
plankton_species <- plankton_data |> 
  filter(!grepl("^Division:|^Class:", taxa, ignore.case = TRUE))


## Extract Species name and abundance data
species_names <- plankton_species$taxa
plankton_data_new <- plankton_species |> 
  select('ibeshe', 'ilashe', 'igboalejo', 'igboeseyero', 'igboologun') |> 
  mutate(across(everything(), ~replace_na(., 0))) # Replace NA with 0

# Convert to matrix and set row names
plankton_matrix <- as.matrix(plankton_data_new)
rownames(plankton_matrix) <- species_names

# Transpose so rows = stations, columns = species
plankton_data_new <- as.data.frame(t(plankton_matrix))

## Create a metadata for the Stations
metadata <- tibble(
  stations = as.factor(c(
    'Ibeshe', 'Ilashe', 'Igboalejo', 'Igboeseyero', 'Igboologun'))
)

## Taxonomic Grouping
taxonomy <- data.frame(
  taxa = species_names,
  Division = NA,
  Class = NA
)

# Fill in division and class info from your original data
current_division <- NA
current_class <- NA


for(i in 1:nrow(plankton_data)) {
  if(grepl("^Division:", plankton_data$taxa[i])) {
    current_division <- gsub("Division: ", "", plankton_data$taxa[i])
  } else if(grepl("^Class:", plankton_data$taxa[i])) {
    current_class <- gsub("Class: ", "", plankton_data$taxa[i])
  } else {
    # This is a species row
    species_idx <- which(taxonomy$taxa == plankton_data$taxa[i])
    if(length(species_idx) > 0) {
      taxonomy$Division[species_idx] <- current_division
      taxonomy$Class[species_idx] <- current_class
    }
  }
}

print(taxonomy)

## Number of Species/Taxa
ncol(plankton_data_new)

## Number of Stations
nrow(plankton_data_new)

## DESCRIPTIVE STATISTICS
total_abundance <- rowSums(plankton_data_new)
print(total_abundance)
taxa_richness <- specnumber(plankton_data_new)
print(taxa_richness)

## Abundance Summary
abundance_summary <- tibble(
  Stations = metadata$stations,
  Total_abundance = total_abundance,
  Number_of_taxa = taxa_richness,
  Mean_Abundance = total_abundance / taxa_richness
)

abundance_summary |> 
  kable()

stats_plankton <- describe(plankton_data_new)

## CALCULATE DIVERSITY INDICES

## Alpha Diversity Indices
diversity_indices <- tibble(
  stations = metadata$stations,
  Shannon = round(diversity(plankton_data_new, index = "shannon"), 2),
  Simpson = round(diversity(plankton_data_new, index = "simpson"), 2),
  InvSimpson = round(diversity(plankton_data_new, index = "invsimpson"), 2),
  Richness = specnumber(plankton_data_new),
  Margalef = round((specnumber(plankton_data_new) - 1) / log(rowSums(plankton_data_new) + 1), 2),
  Pielou_Evenness = round(diversity(plankton_data_new, index = "shannon") / 
    log(specnumber(plankton_data_new)), 2),
  Total_Abundance = rowSums(plankton_data_new)
)

diversity_indices |> 
  kable()


##geom_bar()## Test for Normality of data across all stations
## Using Shapiro Wilk test

shannon_normality <- shapiro.test(diversity_indices$Shannon)
print(shannon_normality)

richness_normality <- shapiro.test(diversity_indices$Richness)
print(richness_normality)

simpson_normality <- shapiro.test(diversity_indices$Simpson)
print(simpson_normality)

margalef_normality <- shapiro.test(diversity_indices$Margalef)
print(margalef_normality)

evenness_normality <- shapiro.test(diversity_indices$Pielou_Evenness)
print(evenness_normality)

abundance_normality <- shapiro.test(diversity_indices$Total_Abundance)
print(abundance_normality)
## Kruskal Wallis Test for Shannon Diversity

shannon_test <- kruskal.test(Shannon ~ stations, data = diversity_indices)

shannon_test 

## Posthoc Test for Shannon Diversity
shannon_posthoc <- dunnTest(Shannon ~ stations, data = diversity_indices)

shannon_posthoc 
## Kruskal Wallis Test for Species Richness

richness_test <- kruskal.test(Richness ~ stations, data = diversity_indices)

richness_test 
## Posthoc Test for Species Richness
richness_posthoc <- dunnTest(Richness ~ stations, data = diversity_indices)

richness_posthoc 

## Kruskal Wallis Test for Simpson Diversity
simpson_test <- kruskal.test(Simpson ~ stations, data = diversity_indices)

simpson_test 

## Posthoc Test for Simpson Diversity
simpson_posthoc <- dunnTest(Simpson ~ stations, data = diversity_indices)

simpson_posthoc 

## Kruskal Wallis Test for Total Abundance
total_ab_test <- kruskal.test(Total_Abundance ~ stations, data = diversity_indices)

total_ab_test

## Posthoc test for Total Abundance
total_posthoc <- dunnTest(Total_Abundance ~ stations, data = diversity_indices)

total_posthoc


### BETA-DIVERSITY ANALYSIS

### PERMANOVA - PERMUTATIONAL MULTIVARIATE ANALYSIS OF VARIANCE 
## test for differences in community composition

### Calculate Bray-Curtis Dissimilarity 
plankton_dis <- vegdist(plankton_data_new, method = 'bray')

plankton_dis |> 
  kable()

## PERMANOVA
set.seed(123)
permanova_result <- adonis2(plankton_data_new ~ metadata$stations,
                            method = 'bray',
                            permutations = 999)

permanova_result |> 
  kable()

## Pair-wise Permanova
pairwise_permanova <- pairwise.adonis(plankton_data_new, 
                                      metadata$stations, 
                                      sim.method = "bray",
                                      p.adjust.m = "bonferroni",
                                      perm = 999)
pairwise_permanova |> 
  kable()


### NMDS ORDINATION
# Perform NMDS for the Sampling Stations

set.seed(123)

nmds_result  <- metaMDS(plankton_data_new, distance = 'bray',
                        trymax = 100, trace = FALSE)

nmds_result

##Extract Site Scores
nmds_score <- as.data.frame(scores(nmds_result, display = 'sites'))
nmds_score$stations <- metadata$stations

station_nmds <- nmds_score |> 
  ggplot(aes(NMDS1, NMDS2, shape = stations, color = stations, label = stations))+
  geom_point(size = 4) +
  geom_text(vjust =1.1, hjust =0.45,  size = 3.5, fontface = "bold") +
  stat_ellipse(aes(group = stations), level = 0.4, linetype = 2, linewidth = 1) +
  theme_bw() +
  labs(title = "NMDS of Phytoplankton Communities",
       subtitle = paste0("Bray-Curtis Similarity\nStress = ", round(nmds_result$stress, 3))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"))
print(station_nmds)

## Perform NMDS Ordination for Species
## (NMDS) ordination of tOP 10 dominant Plankton Diversity
# Extract species scores
species_score <- as.data.frame(scores(nmds_result, display = "species"))
species_score$species <- rownames(species_score)
species_score$abundance <- colSums(plankton_data_new)

# Get top 10 species
top_species <- species_score |> 
  arrange(desc(abundance)) |> 
  head(10)

# NMDS species vectors
species_nmds_10 <- nmds_score |> 
  ggplot(aes(NMDS1, NMDS2))+
  geom_point(aes(color = stations), size = 3) +
  geom_segment(data = top_species, 
               aes(x = 0, y = 0, xend = NMDS1 * 0.8, yend = NMDS2 * 0.8),
               arrow = arrow(length = unit(0.4, "cm")),
               color = "red", alpha = 0.7, linewidth = 0.7) +
  geom_text_repel(data = top_species,
                  aes(x = NMDS1, y = NMDS2, label = species),
                  fontface = "italic", size = 4,
                  box.padding = 0.8, point.padding = 0.5,
                  segment.color = "green", max.overlaps = 20)+
  theme_bw() +
  labs(title = "NMDS with Top 10 Dominant Phytoplankton Diversity",
       subtitle = "Vector length represents species abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
print(species_nmds_10)

# NMDS species vectors

### (NMDS) Ordination of the Sampled Plankton Diversity.
plankton_diversity_nmds <- nmds_score |> 
  ggplot(aes(NMDS1, NMDS2))+
  geom_point(aes(color = stations), size = 3) +
  geom_segment(data = species_score, 
               aes(x = 0, y = 0, xend = NMDS1 * 0.8, yend = NMDS2 * 0.8),
               arrow = arrow(length = unit(0.4, "cm")),
               color = "red", alpha = 0.7, linewidth = 0.7) +
  geom_text_repel(data = species_score,
                  aes(x = NMDS1, y = NMDS2, label = species),
                  fontface = "italic", cex = 3, direction = 'both',
                  box.padding = 0.8, point.padding = 0.5, segment.size = 0.4,
                  segment.color = "green", max.overlaps = 52)+
  theme_bw() +
  labs(title = "Non-Metric Multidimensional(NMDS) Plot of all Phytoplankton Diversity",
       subtitle = "Vector length represents species abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

print(plankton_diversity_nmds)


### CANONICAL CORRESPONDENCE ANALYSIS
# Perform CA
ca_result <- cca(plankton_data_new)

summary(ca_result) 

# Calculate variance explained
ca_var_explained <- ca_result$CA$eig / sum(ca_result$CA$eig) * 100
print(ca_var_explained)

# Extract scores
ca_scores <- as.data.frame(scores(ca_result, display = "sites"))
ca_scores$stations <- metadata$stations

# Plot CA
p_ca <- ggplot(ca_scores, aes(x = CA1, y = CA2, color = stations, shape = stations, label = stations)) +
  geom_point(size = 3) +
  geom_text(vjust =1.1, hjust =0.45,  size = 3.5, fontface = "bold") +
  theme_bw() +
  labs(title = "Correspondence Analysis of Phytoplankton Communities",
       x = paste0("CA1 (", round(ca_var_explained[1], 1), "%)"),
       y = paste0("CA2 (", round(ca_var_explained[2], 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12, face = "bold"))

print(p_ca)


### CLUSTER ANALYSIS

# Hierarchical clustering
plankton_cluster <- hclust(plankton_dis, method = "ave")
plankton_cluS_Comp <- hclust(plankton_dis, method = "complete")


# Plot dendrogram
plot(plankton_cluster, 
     labels = metadata$stations, 
     main = "Hierarchical Clustering of the Phytoplankon Sampling Stations\n(Average Linkage, Bray-Curtis Distance)",
     xlab = "Station", 
     ylab = "Bray-Curtis Dissimilarity",
     sub = "",
     hang = -1,
     cex = 1.2,
     cex.main = 1.3,
     cex.lab = 1.1)
rect.hclust(plankton_cluster, k = 3, border = c("red", "blue", "green"))


## Plot Complete
plot(plankton_cluS_Comp, 
     main = "Hierarchical Clustering of the Phytoplankon Sampling Stations\n(Complete Linkage, Bray-Curtis Distance)",
     xlab = "Station",
     ylab = "Bray-Curtis Dissimilarity",
     hang = -1, cex = 1.2)
rect.hclust(plankton_cluS_Comp, k = 3, border = c("red", "blue", "green"))

# Dissimilarity matrix
dissimilarity_matrix <- as.matrix(plankton_dis)
dissimilarity_matrix <- round(dissimilarity_matrix, 3)
dissimilarity_matrix |> 
  kable()

## BETA DIVERSITY ANALYSIS/PARTITIONING

# Convert to presence-absence
plankton_pa <- decostand(plankton_data_new, method = "pa")

# Partition beta diversity
beta_results <- beta.pair(plankton_pa, index.family = "sorensen")

# Extract components
beta_sim <- as.matrix(beta_results$beta.sim)  # Turnover
beta_sne <- as.matrix(beta_results$beta.sne)  # Nestedness
beta_sor <- as.matrix(beta_results$beta.sor)  # Total

## Turnover Component (Simpson Dissimilarity)
round(beta_sim, 3) |> 
  kable()

## Nestedness Component
round(beta_sne, 3) |> 
  kable()

## Total Beta Diversity (Sorensen Dissimilarity)
round(beta_sor, 3) |> 
  kable()

# Calculate means
mean_turnover <- mean(beta_sim[upper.tri(beta_sim)])
mean_nestedness <- mean(beta_sne[upper.tri(beta_sne)])
mean_total <- mean(beta_sor[upper.tri(beta_sor)])

## Mean Beta Diversity Components:
## Turnover
## Nestedness
## Total

## TAXONOMIC COMPOSITION ANALYSIS

# Transpose data
species_abundance <- as.data.frame(t(plankton_data_new))
species_abundance$taxa <- rownames(species_abundance)

# Join with taxonomy
species_with_taxonomy <- species_abundance |> 
  left_join(taxonomy, by = "taxa")

# Summarize by Phyto-plankton Absolute Abundance according to the Division)
division_abundance <- species_with_taxonomy |> 
  group_by(Division) |> 
  summarise(
    ibeshe = sum(ibeshe, na.rm = TRUE),
    ilashe = sum(ilashe, na.rm = TRUE),
    igboalejo = sum(igboalejo, na.rm = TRUE),
    igboeseyero = sum(igboeseyero, na.rm = TRUE),
    igboologun = sum(igboologun, na.rm = TRUE),
    .groups = "drop"
  )

## Absolute Abundance by Division
division_abundance |> 
  kable()

# Calculate relative abundance of Phyto-plankton Communities by Division
division_relative <- division_abundance
station_cols <- c("ibeshe", "ilashe", 
                  "igboalejo", "igboeseyero", "igboologun")

division_relative <- division_relative |> 
  mutate(
    ibeshe = (ibeshe / sum(ibeshe)) * 100,
    ilashe  = (ilashe  / sum(ilashe )) * 100,
    igboalejo = (igboalejo / sum(igboalejo)) * 100,
    igboeseyero = (igboeseyero / sum(igboeseyero)) * 100,
    igboologun = (igboologun / sum(igboologun)) * 100
  )

## Relative Abundance by Division
division_relative |> 
  kable()


# Visualize Division composition
division_plot_data <- division_relative |> 
  pivot_longer(cols = all_of(station_cols), 
               names_to = "Station", 
               values_to = "Relative_Abundance")

p_division <- division_plot_data |> 
  ggplot(aes(x = Station, y = Relative_Abundance, fill = Division)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black",
           linewidth = 0.5, width = 0.7) +
  theme_bw() +
  labs(title = "Phytoplankton Communities Composition by Division",
       y = "Relative Abundance (%)",
       x = "Station") +
  scale_fill_brewer(palette = "Set2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13, face = "bold"),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11))

print(p_division)

## Summarize the Absolute Abundance of Phytoplankton Communities by Class
class_abundance <- species_with_taxonomy |> 
  group_by(Division, Class) |> 
  summarise(
    ibeshe = sum(ibeshe, na.rm = TRUE),
    ilashe = sum(ilashe, na.rm = TRUE),
    igboalejo = sum(igboalejo, na.rm = TRUE),
    igboeseyero = sum(igboeseyero, na.rm = TRUE),
    igboologun = sum(igboologun, na.rm = TRUE),
    .groups = "drop"
  ) 

class_abundance |> 
  kable()

## Calculate relative abundance of Phytoplankton Communities by class
class_relative <- class_abundance |> 
  mutate(
    ibeshe = (ibeshe / sum(ibeshe)) * 100,
    ilashe  = (ilashe  / sum(ilashe )) * 100,
    igboalejo = (igboalejo / sum(igboalejo)) * 100,
    igboeseyero = (igboeseyero / sum(igboeseyero)) * 100,
    igboologun = (igboologun / sum(igboologun)) * 100
  )

class_relative |> 
  kable()

## Visualize Class composition
class_plot_data <- class_relative |> 
  pivot_longer(cols = all_of(station_cols), 
               names_to = "Station", 
               values_to = "Relative_Abundance")

p_class <- ggplot(class_plot_data, 
                  aes(x = Station, y = Relative_Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black",
           linewidth = 0.5, width = 0.7) +
  theme_bw() +
  labs(title = "Phytoplankton Community Composition by Class",
       y = "Relative Abundance (%)",
       x = "Station") +
  scale_fill_manual(values = rainbow(length(unique(class_relative$Class)))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.text = element_text(size = 9))

print(p_class)

# Count taxa per division
taxa_per_division <- taxonomy |> 
  group_by(Division) |> 
  summarise(
    Number_of_Taxa = n(),
    Number_of_Classes = n_distinct(Class)
  ) |> 
  arrange(desc(Number_of_Taxa))

taxa_per_division |> 
  kable()

# Count taxa per class
taxa_per_class <- taxonomy |> 
  group_by(Division, Class) |> 
  summarise(Number_of_Taxa = n(), .groups = "drop") |> 
  arrange(Division, desc(Number_of_Taxa))

taxa_per_class |> 
  kable()

## Environmental Data

## Import the data set for Analysis
water_quality <- read_excel('Water Quality Analysis(Mr. Bosun).xlsx')

## Data Preparation
water_quality <- water_quality |> 
  janitor::clean_names() |> 
  mutate(
    stations = as.factor(stations),
    season = as.factor(season),
    months = as.factor(months)
  )

water_quality$months <- factor(water_quality$months, 
                               levels = unique(water_quality$months))
water_quality$stations <- factor(water_quality$stations,
                                 levels = unique(water_quality$stations))

## Mean of Physicochemical parameters
water_mean <- water_quality |> 
  select(stations, dissolved_oxygen, p_h, salinity, temperature,
         total_dissolved_solid, turbidity) |> 
  group_by(stations) |> 
  summarise(
    Dissolved_Oxygen = round(mean(dissolved_oxygen), 2),
    pH = round(mean(p_h), 2),
    Salinity = round(mean(salinity), 2),
    Temperature = round(mean(temperature), 2),
    Total_Dissolved_Solids = round(mean(total_dissolved_solid), 2),
    Turbidity = round(mean(turbidity), 2),
    .groups = 'drop'
  )


# ==========================================
# STEP 2: PERFORM CANONICAL CORRESPONDENCE ANALYSIS
# ==========================================
# Run CCA
cca_result <- cca(plankton_data_new ~ Temperature
                  + Salinity + pH + Dissolved_Oxygen + 
                    Total_Dissolved_Solids + Turbidity,
                  data = water_mean)

# Print summary

print(summary(cca_result))

# Test significance

# Overall test
set.seed(123)
cca_anova <- anova(cca_result, permutations = 999)
print(cca_anova)

# Test by axis
cca_axis <- anova(cca_result, by = "axis", permutations = 999)
print(cca_axis)

# Test by term: Significance by environmental variable
cca_terms <- anova(cca_result, by = "terms", permutations = 999)
print(cca_terms)

# Variance explained
variance_explained <- eigenvals(cca_result) / sum(eigenvals(cca_result)) * 100
print(variance_explained[1:2])

`# ==========================================
# STEP 3: EXTRACT SCORES FOR PLOTTING
# ==========================================

# Extract site scores (stations)
site_scores <- as.data.frame(scores(cca_result, display = "sites", choices = c(1, 2)))
site_scores$Station <- rownames(site_scores)
site_scores$Type <- "Site"

# Extract species scores (phytoplankton)
species_scores <- as.data.frame(scores(cca_result, display = "species", choices = c(1, 2)))
species_scores$Species <- rownames(species_scores)
species_scores$Type <- "Species"
species_scores$Abundance <- colSums(plankton_data_new)

plankton_data_new <- plankton_data_new |> 
  select(-"Tabellaria fenestrata")


# Extract environmental variable scores (arrows)
env_scores <- as.data.frame(scores(cca_result, display = "bp", choices = c(1, 2)))
env_scores$Variable <- rownames(env_scores)
env_scores$Type <- "Environment"

# Get top 10 species for clarity
top_species <- species_scores %>%
  arrange(desc(Abundance)) %>%
  head(10) 

# Add taxonomy to species
species_scores <- species_scores |> 
  left_join(taxonomy, by = c("Species" = "taxa"))

top_species <- top_species |> 
  left_join(taxonomy, by = c("Species" = "taxa"))


## TRIPLOT WITH SPECIES AND environmental variables
# Station colors
station_colors <- c("ibeshe" = "#4DBBD5", 
                    "ilashe" = "#E64B35",
                    "igboalejo" = "#00A087",
                    "igboeseyero" = "#F39B7F",
                    "igboologun" = "#8491B4")
# Scaling factor for environmental arrows
arrow_scale <- 1.5

# Create triplot showing top 10 dominant species
p_triplot <- ggplot() +
  # Add species points (background, gray)
  geom_point(data = species_scores,
             aes(x = CCA1, y = CCA2),
             size = 2.0, color = "green", alpha = 0.3) +
  
  # Add site points (stations)
  geom_point(data = site_scores,
             aes(x = CCA1, y = CCA2, fill = Station),
             size = 8, shape = 21, color = "black", stroke = 1.5) +
  # Add site labels
  geom_text(data = site_scores,
            aes(x = CCA1, y = CCA2, label = Station),
            fontface = "bold", size = 4.5, vjust = -1.0) +
  # Add environmental arrows
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, 
                   xend = CCA1 * arrow_scale, 
                   yend = CCA2 * arrow_scale),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "red3", size = 0.5, alpha = 0.8) +
  # Add environmental variable labels
  geom_text_repel(data = env_scores,
                  aes(x = CCA1 * arrow_scale, 
                      y = CCA2 * arrow_scale, 
                      label = Variable),
                  color = "red3", fontface = "bold", size = 3,
                  box.padding = 0.5, point.padding = 0.3,
                  segment.color = "red3", segment.size = 0.3) +
  # Add top 10 species labels
  geom_text_repel(data = top_species,
                  aes(x = CCA1, y = CCA2, label = Species),
                  color = "black", fontface = "italic", size = 3,
                  box.padding = 0.3, point.padding = 0.2,
                  segment.color = "gray50", segment.size = 0.2,
                  max.overlaps = 15) +
  # Scales and theme
  scale_fill_manual(values = station_colors) +
  theme_bw() +
  labs(title = "CCA Triplot: Phytoplankton-Environment Relationships",
       subtitle = paste0('Canonical Correspondence Analysis Showing the Distribution
                         of Top 10 Species' ),
       x = paste0("CCA1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("CCA2 (", round(variance_explained[2], 1), "%)"),
       fill = "Station") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
print(p_triplot)


## TRIPLOT WITH SPECIES BY DIVISION

p_triplot_taxonomy <- ggplot() +
  # Add species points colored by Division
  geom_point(data = species_scores,
             aes(x = CCA1, y = CCA2, colour = Division),
             size = 2, alpha = 0.6) +
  # Add site points
  geom_point(data = site_scores,
             aes(x = CCA1, y = CCA2, fill = Station),
             size = 4, shape = 21, color = "black", stroke = 2) +
  # Add site labels
  geom_text(data = site_scores,
            aes(x = CCA1, y = CCA2, label = Station),
            fontface = "bold", size = 3, vjust = -2.2) +
  geom_text_repel(data = species_scores,
                  aes(x = CCA1, y = CCA2, label = Species),
                  fontface = "italic", size = 3.5,
                  box.padding = 0.4, max.overlaps = 60) +
  # Add environmental arrows
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, 
                   xend = CCA1 * arrow_scale, 
                   yend = CCA2 * arrow_scale),
               arrow = arrow(length = unit(0.5, "cm")),
               color = "blue", size = 1.0, linewidth = 0.7) +
  # Add environmental labels
  geom_text_repel(data = env_scores,
                  aes(x = CCA1 * arrow_scale, 
                      y = CCA2 * arrow_scale, 
                      label = Variable),
                  color = "red3", fontface = "bold", size = 3,
                  box.padding = 0.5) +
  # Scales
  scale_fill_manual(values = station_colors) +
  scale_color_brewer(palette = "Set1", name = "Phytoplankton\nDivision") +
  theme_bw() +
  labs(title = "Canonical Correspondence Analysis Triplot: Species Colored by Taxonomic Division",
       subtitle = paste0("Showing relationship between phytoplankton communities and environmental variables"),
       x = paste0("CCA1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("CCA2 (", round(variance_explained[2], 1), "%)"),
       fill = "Station") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "right")

print(p_triplot_taxonomy)

# SEPARATE BIPLOTS (SITES + ENV, SPECIES + ENV)
# Biplot 1: Sites and Environment
p_biplot_sites <- ggplot() +
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, xend = CCA1 * arrow_scale, yend = CCA2 * arrow_scale),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "red3", size = 1.2, linewidth = 0.7) +
  geom_point(data = site_scores,
             aes(x = CCA1, y = CCA2, fill = Station),
             size = 8, shape = 21, color = "black", stroke = 2) +
  geom_text(data = site_scores,
            aes(x = CCA1, y = CCA2, label = Station),
            fontface = "bold", size = 4) +
  geom_text_repel(data = env_scores,
                  aes(x = CCA1 * arrow_scale, y = CCA2 * arrow_scale, label = Variable),
                  color = "red3", fontface = "bold", size = 4.5,
                  box.padding = 0.5) +
  scale_fill_manual(values = station_colors) +
  theme_bw() +
  labs(title = "CCA Biplot: Stations and Environmental Variables",
       x = paste0("CCA1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("CCA2 (", round(variance_explained[2], 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right")

print(p_biplot_sites)

# Biplot 2: Species and Environment
p_biplot_species <- ggplot() +
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, xend = CCA1 * arrow_scale, yend = CCA2 * arrow_scale),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "red3", size = 1.2, linewidth = 0.5) +
  geom_point(data = species_scores,
             aes(x = CCA1, y = CCA2, color = Division),
             size = 2.5, alpha = 0.6) +
  geom_text_repel(data = species_scores,
                  aes(x = CCA1, y = CCA2, label = Species),
                  fontface = "italic", size = 3.5,
                  box.padding = 0.4, max.overlaps = 15) +
  geom_text_repel(data = env_scores,
                  aes(x = CCA1 * arrow_scale, y = CCA2 * arrow_scale, label = Variable),
                  color = "red3", fontface = "bold", size = 4,
                  box.padding = 0.5) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "CCA Biplot: Species and Environmental Variables",
       subtitle = "Phytoplankton species labeled",
       x = paste0("CCA1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("CCA2 (", round(variance_explained[2], 1), "%)"),
       color = "Division") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "right")

print(p_biplot_species)