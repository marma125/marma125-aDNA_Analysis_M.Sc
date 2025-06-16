# Description ------------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Author: Max Marklein
# University of Bremen
# Last updated: 16.06.2025
# 
# The script (Msc_Analysis.R) curates the data and performs the analysis discussed
# in the manuscript. It builds upon previous parsing scripts written by
# Dr. Chiara Vanni, therefore the original data sources where not manipulated.
# 
# Data files used:
# 1. hires-organelles-viruses-smags.tax.tsv
# 2. sample_table_cdata.tsv
# 3. pilot_cores_aDNA.txt
# 4. List_Samples_aDNA-Atlantic_with_age.xlsx
# 5. tp-mapping-filtered.nocontam.all_samples.smpl_agg.damage_filt_age.tax.tsv.gz
# 6. tp-mdmg.weight-1.local.all_samples.smpl_agg.damage_filt_age_euk_cyano_arc.csv.gz
# 7. aDNA_APCI_Proxy.csv
# 8. aDNA_TOC.csvincluding all TSS scaled, damaged and mapped reads from all samples.
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Set Up -----------------------------------------------------------------------

## Set Working directory ----
setwd("~/Wichtiges/Studium/Thesis M.Sc/marma125-aDNA_Analysis_M.Sc")

## clear environment ----
rm(list=ls())
dev.new()

## Load necessary libraries ----
library(tidyverse)
library(patchwork)
library(ggrepel)
library(shadowtext) 
library(data.table)
library(vegan)
#library(RColorBrewer)
library(paletteer)
#library(factoextra)
library(WGCNA)
library(pheatmap)
library(igraph)
library(mixOmics)
library(treemap)
library(readxl)
library("rnaturalearth")
library("rnaturalearthdata")

# Coustom colors for coloring the cores
colorful <- c(
  "GeoB18549-2" = "#E69F00FF",  # orange
  "GeoB16320-2" = "#56B4E9FF",  # sky blue
  "GeoB3938-1" = "#009E73FF",  # bluish green
  "GeoB1706-2" = "#BD1719",  # yellow
  "GeoB6408-4" = "#0072B2FF",  # blue
  "GeoB1903-3" = "#D55E00FF",  # vermilion
  "GeoB2004-2" = "#CC79A7FF",  # reddish purple
  "GeoB9506-1" = "#999999FF",  # gray
  "GeoB2215-10" = "#A6761DFF",  # brown
  "GeoB9064-1" = "#F0E442FF"   # teal green
)

#______________________________________________________________________________#

# Load data --------------------------------------------------------------------

## Taxa data ----

tax_info <- read_tsv("taxonomy_metadata/hires-organelles-viruses-smags.tax.tsv",
                     show_col_types = FALSE,
                     col_names = c("reference", "tax_string")
) %>%
  separate( 
    col = tax_string,
    sep = ";",
    into = c(
      "domain",
      "lineage",
      "kingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species",
      "strain"
    )
  )

# Sample contextual data
samples <- read_tsv("sample_table_cdata.tsv", show_col_types = FALSE)

# Core and age list
core_list <- unique(str_subset(samples$core, "^GeoB"))
ages <- sort(unique(samples$age))

# Core geographic information
cores_map <- fread("core_metadata/pilot_cores_aDNA.txt") %>%
  mutate(latitude=as.numeric(latitude)) %>%
  mutate(year=as.character(year)) %>%
  inner_join(samples %>% dplyr::select(core) %>% distinct()) %>%
  mutate(core_labels=paste0(core,"; ",year)) %>%
  mutate(core_labels=fct_reorder(core_labels,as.numeric(year))) %>% 
  mutate(core = fct_reorder(as.factor(core), abs(latitude)))

# Age model
real_ages <- read_excel("core_metadata/List_Samples_aDNA-Atlantic_with_age.xlsx") %>% 
  dplyr::select("GeoB core", Depth, Age) %>% 
  rename(core = "GeoB core", r_age = Age) %>% 
  drop_na() %>% 
  distinct() %>% 
  group_by(core) %>% 
  mutate(age = case_when(
    core == "GeoB1706-2"        ~ seq(14, 14 + 2 * (n() - 1), by = 2),
    n() == 9                    ~ seq(10, 10 + 2 * (n() - 1), by = 2),
    TRUE                        ~ seq(8, 8 + 2 * (n() - 1), by = 2)
  )) %>%
  ungroup()

# Taxa data aggregated at the sample level
tax_data_agg <- read_tsv("results/tp-mapping-filtered.nocontam.all_samples.smpl_agg.damage_filt_age.tax.tsv.gz", show_col_types = FALSE)

# Damage data aggregated at sample level
dmg_local_agg <- read_csv("results/tp-mdmg.weight-1.local.all_samples.smpl_agg.damage_filt_age_euk_cyano_arc.csv.gz", show_col_types = FALSE)

# taxonomic annotations with damage
tax_data_agg <- tax_data_agg %>%
  inner_join(dmg_local_agg) %>%
  mutate(label=paste(core_labels,age,sep=":"))

# Get relative abundances by Total sum scaling (TSS) and filter out Eukaryota and non damaged taxa
tax_rel_abund_genus <- tax_data_agg %>%
  inner_join(tax_info) %>%
  filter(domain %in% c("d__Archaea","d__Bacteria")) %>%
  #filter(domain %in% c("d__Eukaryota")) %>%
  #filter(core != "GeoB1706-2") %>% # Remove Namibian core - Upwelling environment
  dplyr::select(label, core_labels,age, domain, phylum, genus, abundance, is_dmg, reference) %>%
  group_by(label, core_labels,age, domain, phylum, genus, is_dmg) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  group_by(label,core_labels,age) %>%
  filter(is_dmg == "Damaged") %>% #keeping both bacteria and archaea that are damaged
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(
    abundance = abundance / total_abundance
  ) %>%
  mutate(genus = ifelse(abundance < 0.1, "Other", genus)) %>%
  group_by(label, core_labels,age, domain, genus, is_dmg) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(genus=gsub("g__","",genus)) %>%
  mutate(genus = fct_reorder(genus, abundance, "sum", .desc = TRUE)) %>%
  mutate(genus=case_when(domain=="d__Bacteria" ~ paste("B",genus,sep="_"),
                         domain=="d__Archaea" ~ paste("A",genus,sep="_"),
                         domain=="d__Eukaryota" ~ paste("E",genus,sep="_"),
                         TRUE ~ paste("V",genus,sep="_"))) %>%
  inner_join(cores_map %>% dplyr::select(core_labels, latitude, longitude)) %>%
  separate(core_labels, into = c("core", "year"), sep = ";") %>%
  mutate(core = fct_reorder(as.factor(core), abs(latitude))) %>% 
  inner_join(real_ages) %>% 
  dplyr::select(-age, -Depth) %>% 
  rename(age = r_age)


# Make proportions to 1 and convert to wide format matrix with colnames having the sample id
otu_tss <- tax_rel_abund_genus %>% dplyr::select(-domain) %>% distinct() %>%
  mutate(abundance=abundance) %>%
  ungroup() %>% distinct() %>% 
  pivot_wider(names_from = "genus", values_from = "abundance", values_fill = 0) %>% 
  column_to_rownames("label") # Convert to matrix format (first column have sample ID)

sum(otu_tss[1, -c(1:6)])

## Environmental data ----

cdata <- read_delim("core_metadata/aDNA_APCI_Proxy.csv",delim = ";", show_col_types = F)

EA <- read_delim("core_metadata/aDNA_TOC.csv",delim = ";", show_col_types = F) %>% 
  rename(depth = Depth_cm, core = Core_ID)

cdata_all <- cdata %>% 
  dplyr::select(Core_ID,Depth_cm,`UK37_SST_°C`,
                `TEXh86_0-200_°C`,`CCaT_SST_°C`, year) %>%
  mutate( `UK37_SST_°C`=as.numeric(`UK37_SST_°C`),
          `TEXh86_0-200_°C` = as.numeric(`TEXh86_0-200_°C`),
          `CCaT_SST_°C` = as.numeric(`CCaT_SST_°C`)) %>%
  rename(depth = Depth_cm, core = Core_ID) %>% 
  inner_join(samples, by = c("core","depth")) %>% 
  dplyr::select(core, year, age, depth, `UK37_SST_°C`) %>% 
  distinct() %>% 
  inner_join(EA, by = c("core","depth")) %>% 
  mutate(core_labels=paste0(core,"; ",year)) %>%
  mutate(core_labels=fct_reorder(core_labels,as.numeric(year))) %>% 
  mutate(label=paste(core_labels,age,sep=":")) %>% 
  dplyr::select(where(~ all(!is.na(.)))) %>% 
  column_to_rownames("label") %>% 
  dplyr::select(core, year, age, `UK37_SST_°C`, TOC) %>% 
  inner_join(real_ages) %>% 
  dplyr::select(-age, -Depth) %>% 
  rename(age = r_age)

## Combine OTU and Env data
otu_tss_env <- inner_join(otu_tss, cdata_all, by = c("core", "age")) %>% 
  dplyr::select(-year.x, -year.y, -is_dmg) %>%
  dplyr::select(c(core, age, latitude, longitude, `UK37_SST_°C`, TOC), everything()) %>% 
  mutate(label = paste(core, age, sep = "_")) %>% 
  column_to_rownames("label")

#______________________________________________________________________________#

# Fig 1: Data Introduction ----------------------------------------------------
## Plotting loaded data ----
### A: Map ----

a <- ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(fill="#272733") +
  theme_bw() +
  coord_sf(xlim = c(-95, 35), expand = FALSE) +
  geom_jitter(data = cores_map,
              aes(x = longitude, y = latitude, fill = core),
              size=5.5, shape=21, alpha=0.8, position="jitter") +
  scale_fill_manual(values = colorful) +
  #paletteer::scale_fill_paletteer_d("ggthemr::copper") +
  geom_label_repel(data=cores_map,
                   aes(x=longitude, 
                       y=latitude,
                       label=core),
                   color="black", size = 2.5, fill = "white", box.padding = 0.39) +
  labs(x = "Longitude", y = "Latitude", fill = "Core") +
  theme(legend.position = "none")


### B: Env data ----

b <- ggplot(data = cdata_all, aes(x = age)) +
  geom_point(aes(y = `UK37_SST_°C`, group = core), color = "#CFB8B9", alpha = 0.7) +
  geom_line(aes(y = `UK37_SST_°C`, group = core), color = "#CFB8B9", alpha = 0.5) +
  geom_point(aes(y = TOC * 5, group = core), color = "#C1C0DB", alpha = 0.7) +
  geom_line(aes(y = TOC * 5, group = core), color = "#C1C0DB", alpha = 0.5) +
  geom_smooth(aes(y = `UK37_SST_°C`, color = "UK'37 SST [°C]"), method = "loess", se = F) +
  geom_smooth(aes(y = TOC * 5, color = "TOC [%]"), method = "loess", se = F) +
  scale_x_continuous(breaks=seq(2,30,by=4)) +
  scale_y_continuous(
    name = "UK'37 SST [°C]",
    sec.axis = sec_axis(~./5, name = "TOC [%]")) +
  scale_color_manual(name = "Loess Fit", 
                     values = c("UK'37 SST [°C]" = "#C8102EFF", "TOC [%]" = "#0C2340FF")) +
  theme_bw() +
  labs(x = "Age [ka BP]") +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA))

### C: Relative Abundance ----

c <- ggplot(tax_rel_abund_genus, aes(x = age, y = abundance, fill = genus)) +
  geom_col(position = "stack", color = "black", linewidth = 0.3, width = 1) +
  theme_bw() +
  xlab("Age [ka BP]") +
  ylab("Relative Abundance [%]") +
  facet_wrap(~core) +
  #scale_fill_manual(name="genus", values=paletteTax(length(unique(tax_rel_abund_genus$genus)))) +
  paletteer::scale_fill_paletteer_d("ggsci::default_igv", direction = 1) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks=seq(2,30,by=4)) +
  labs(fill = "Genera") +
  theme_bw() +
  guides(
    fill = guide_legend(
      keywidth = 0.6,     # Width of legend rectangles (cm)
      keyheight = 0.4,    # Height of legend rectangles (cm)
      default.unit = "cm"
    )
  ) +
  theme(
    legend.position = "top",
    legend.key = element_rect(size = 0.2),              # Smaller key boxes
    legend.key.height = unit(0.3, "cm"),                # Reduce vertical size
    legend.key.width = unit(0.5, "cm"),                 # Reduce horizontal size
    legend.text = element_text(size = 9),               # Smaller text
    legend.title = element_text(size = 12),              # Smaller legend title
    legend.spacing.y = unit(0.1, 'cm'),                 # Reduce vertical gap
    legend.margin = margin(0, 0, 0, 0),                 # Remove outer margin
    legend.box.margin = margin(0, 0, 0, 0)              # Remove box margin
  )

## Patch ----

(a + b) / c + plot_annotation(tag_levels = "A")

ggsave("Plots/Fig1.pdf", width = 27.75, height = 30, units = "cm")

# Fig 2: Biodiversity ----------------------------------------------------------

## Alpha Diversity ----
# make proportions to 1 and covert to wide format
otu_tss_smpl <- tax_rel_abund_genus %>% dplyr::select(-domain) %>% distinct() %>%
  ungroup() %>% distinct() %>%
  pivot_wider(names_from = "genus", values_from = "abundance", values_fill = 0)

otu_tss_smpl_wide <- tax_rel_abund_genus %>% dplyr::select(-domain) %>% distinct() %>%
  ungroup() %>% distinct() %>%
  pivot_wider(names_from = "genus", values_from = "abundance", values_fill = 0) %>%
  arrange(label) %>%
  dplyr::select(-c(2:7)) %>%
  as.data.frame(.) %>%
  column_to_rownames("label")

# Test the significance of differences in composition among years. 
# Interpretation: The changes in composition over time were significant and different transects changed differently
# using the bray curties dissimilarity (Beta diversity)
adonis2(otu_tss_smpl_wide ~ otu_tss_smpl$core * otu_tss_smpl$age, method = 'bray')

comm_diversity <- function(id, data) {
  
  # Filter by core
  df_work <- dplyr::filter(data, core == id)
  
  # Subset just the species (assumes these start at col 8)
  species_data <- df_work[-c(1:7)]
  
  # Calculate diversity indices
  H <- vegan::diversity(species_data) # Shannon
  specnum <- vegan::specnumber(species_data) # Richness
  S <- vegan::diversity(species_data, index = "invsimpson") # Inverse Simpson
  
  # Add indices to table
  df_work <- df_work %>%
    mutate(Shannon = H, 
           iSimpson = S,
           Richness = specnum)
  
  # Richness change per year (slope of Richness ~ age)
  richness_lm <- lm(Richness ~ age, data = df_work)
  richness_slope <- coef(richness_lm)["age"]
  
  # Species gains/losses (relative to oldest sample)
  ref_sample <- df_work %>% filter(age == max(age)) %>% dplyr::select(-(1:7)) %>% colSums()
  later_samples <- df_work %>% filter(age != max(age)) %>% dplyr::select(-(1:7)) %>% colSums()
  
  species_gained <- sum((later_samples > 0) & (ref_sample == 0))
  species_lost <- sum((later_samples == 0) & (ref_sample > 0))
  
  df_work <- df_work %>%
    mutate(Richness_slope = richness_slope,
           Species_gains = species_gained,
           Species_losses = species_lost)
  
  return(df_work)
}

# Recalculate the Alpha diversity measures
comm_div <- lapply(core_list, comm_diversity, data = otu_tss_smpl) %>% 
  bind_rows() %>% 
  dplyr::select(label, core, age, latitude, longitude, Shannon, iSimpson, Richness,
                Richness_slope, Species_gains, Species_losses) %>% 
  inner_join(otu_tss_env, by = c("core", "age", "latitude", "longitude")) %>% 
  column_to_rownames("label") %>% 
  group_by(core) %>% 
  mutate(Richness_slope_average = mean(Richness_slope)) %>% 
  ungroup() %>% 
  dplyr::select(c(1:10), "Richness_slope_average")

## Beta Diversity (Dispersion) ----
# Are the individual centroids of the core in the BC distance matrix significantly different?
# Calculate dissimilarity matrix (Bray Curtis) (Beta diversity measure)
bc_dist <- vegdist(otu_tss[-(1:6)], method = "bray")
dispersion <- betadisper(bc_dist, group=otu_tss$core, type = "centroid")
anova(dispersion)
#plot(dispersion)
adonis2(bc_dist~as.factor(otu_tss$core), data=bc_dist, permutations=9999)
dispersion.res <- as.data.frame(scores(dispersion)$sites) %>%
  tibble::rownames_to_column("label") %>%
  mutate(core = dispersion$group) %>%
  left_join(otu_tss_smpl[c(1:7)]) %>% 
  dplyr::select(-is_dmg, -year) %>% 
  mutate(core = fct_reorder(as.factor(core), abs(latitude)))

PCoA.hulls <- dispersion.res %>%
  group_by(core) %>%
  slice(chull(PCoA1, PCoA2)) %>%  # Get hull vertices
  ungroup() %>% 
  mutate(core = fct_reorder(as.factor(core), abs(latitude)))

# 64.4% of the variance in community composition (based on Bray-Curtis distances) 
# is explained by the grouping factor core.
# F = 17.122 and p = 1e-04: The difference between the groups is highly significant, 
# meaning there are clear composition differences between the core groups

### A: Richness Map ----
a <- ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(fill="#272733") +
  theme_bw() +
  coord_sf(xlim = c(-95, 35), expand = FALSE) +
  geom_point(data = comm_div,
              aes(x = longitude, y = latitude, fill = Richness_slope_average),
              size=5.5, shape=21, alpha=0.8) +
  scale_fill_gradient2(midpoint = 0, low = "darkblue", mid = "white", high = "darkred") +
  # geom_label_repel(data=cores_map,
  #                  aes(x=longitude, 
  #                      y=latitude,
  #                      label=core),
  #                  color="black", size = 2.5, fill = "white", box.padding = 0.39) +
  labs(x = "Longitude", y = "Latitude", fill = "Average richness change\n(no. of species per ka)")

### B: Density ----
b <- ggplot(comm_div, aes(x = Richness_slope)) +
  geom_density(colour = "black", fill = "black", alpha = 0.3) +
  geom_vline(aes(xintercept = mean(Richness_slope, na.rm = TRUE)),
             colour = "black", linetype = "dashed", linewidth = 0.5) +
  # geom_vline(aes(xintercept = 0),
  #            colour = "grey40", linetype = "dotted", size = 0.5) +
  labs(x = "Richness change", y = "Density") +
  theme_bw()

### C: PCoA Metric Dispersion analysis and Multidimensional scaling ----
c <-
  ggplot() +
  geom_polygon(data = PCoA.hulls, aes(x = PCoA1, y = PCoA2, fill = core, color = core), linetype = "dashed", alpha = 0.2) +
  geom_point(data = dispersion.res, aes(x = PCoA1, y = PCoA2, fill = core, size = age), shape = 21, color = "black") +
  geom_shadowtext(data = dispersion.res %>% 
              group_by(core) %>% 
              summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2)),
            aes(x = PCoA1, y = PCoA2, label = core), size = 3, color = "black", bg.color = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = scales::alpha("black", 0.15)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = scales::alpha("black", 0.15)) + 
  scale_color_manual(values = colorful, name = "Core") +
  scale_fill_manual(values = colorful, name = "Core") +
  #paletteer::scale_fill_paletteer_d("ggthemr::copper", name = "Core") +
  #paletteer::scale_color_paletteer_d("ggthemr::copper", name = "Core") +
  labs(
    x = paste0("PCoA1 (", round(dispersion$eig[1] / sum(dispersion$eig) * 100, 1), "%)"),
    y = paste0("PCoA1 (", round(dispersion$eig[2] / sum(dispersion$eig) * 100, 1), "%)"),
    size = "Age [ka BP]") +
  theme_bw() +
  theme(legend.position = "right")

## Patch ----

(b + a) / c + plot_annotation(tag_levels = "A")

ggsave("Plots/Fig2.pdf", width = 27.75, height = 30, units = "cm")

# Fig 3: CCA -------------------------------------------------------------------
head(otu_tss_env[-(1:6)]) # abundance data
head(otu_tss_env[c(2,5,6)]) # environmental data: Age, SST, TOC only
env_scaled <- data.frame(scale(otu_tss_env[c(2,5,6)])) # Z score normalization: mean of 0 and standard deviation of 1

#CCA on full dataset
res.cca <- cca(otu_tss_env[-(1:6)] ~ ., # abundance data
               data = env_scaled) # environmental data

summary(res.cca)

a <- anova(res.cca) # analysis of variance
a

cca.var_expl <- summary(res.cca)$constr.chi / res.cca$tot.chi * 100
axis_labels <- paste0("CCA", 1:2, " (", round(cca.var_expl[1:2], 1), "%)")

anova(res.cca, by = "axis")
# the first 4 CCA axis are significant

a <- anova(res.cca, by = "term")
a

100 - a["Residual","ChiSquare"] / sum(a$ChiSquare) * 100

anova(res.cca, by = "margin")

# Get scores from the CCA object
cca.site_scores <- as.data.frame(vegan::scores(res.cca, display = "sites"))
cca.site_scores$label <- rownames(cca.site_scores)
cca.site_scores <- cca.site_scores %>% 
  inner_join(rownames_to_column(otu_tss_env, var = "label"), join_by(label)) %>% 
  dplyr::select(c(1:9))
cca.species_scores <- as.data.frame(vegan::scores(res.cca, display = "species"))
cca.biplot_scores <- as.data.frame(vegan::scores(res.cca, display = "bp"))

cca.hulls <- cca.site_scores %>%
  group_by(core) %>%
  slice(chull(CCA1, CCA2)) %>%
  ungroup()

#Plot
ggplot() +
  geom_polygon(data = cca.hulls, aes(x = CCA1, y = CCA2, fill = core, color = core), linetype = "dashed", alpha = 0.2) +
  geom_point(data = cca.site_scores, aes(x = CCA1, y = CCA2, size = age, fill = core), shape = 21, color = "black") +
  geom_shadowtext(data = cca.site_scores %>% 
                    group_by(core) %>% 
                    summarise(CCA1 = mean(CCA1), CCA2 = mean(CCA2)),
                  aes(x = CCA1, y = CCA2, label = core), size = 3, color = "black", bg.color = "white") +
  geom_segment(data = cca.biplot_scores, aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkblue") +
  geom_text(data = cca.biplot_scores, aes(x = CCA1, y = CCA2, label = rownames(cca.biplot_scores)),
            color = "darkblue", size = 4, vjust = 1.3, hjust = -0.1) +
  geom_point(data = cca.species_scores, aes(x = CCA1, y = CCA2), color = "grey", shape = 3, size = 1.2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = scales::alpha("black", 0.15)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = scales::alpha("black", 0.15)) + 
  scale_color_manual(values = colorful, name = "Core") +
  scale_fill_manual(values = colorful, name = "Core") +
  labs(
    x = paste0("CCA1 (", round(res.cca$CCA$eig[1] / res.cca$tot.chi * 100, 1), "%)"),
    y = paste0("CCA2 (", round(res.cca$CCA$eig[2] / res.cca$tot.chi * 100, 1), "%)"),
    size = "Age [ka BP]") +
  theme_bw() +
  theme(legend.position = "right")

ggsave("Plots/Fig3.pdf", width = 27.75, height = 20, units = "cm")

# Fig 4: sPLS ------------------------------------------------------------------
env_scaled <- data.frame(scale(otu_tss_env[c(2:6)])) # include lat/ long


res.sPLS <- spls(X = otu_tss_env[-(1:6)],
                 Y = env_scaled,
                 ncomp = 2, # two PC axis
                 mode = "regression")

explained_variance(res.sPLS$Y, res.sPLS$variates$Y, ncomp = 2)
# comp1     comp2 
# 0.3903787 0.1996191 
# All together I can explain 59% of the variance with my environmental parameters

pdf("Plots/Fig4.pdf", width = 10, height = 13)
cim(res.sPLS, 
    comp = c(1,2), # visualizes the combined loading of the first two axis
    margins = c(10, 13),
    #color.legend = "Correlation",
    row.names = F,
    #title = "59% var. explained by the first 2 components"
)

dev.off()

# Fig 5: WGCNA Summary ---------------------------------------------------------
## Network creation ------------------------------------------------------------
head(otu_tss_env[-(1:6)]) # abundance data TSS
head(otu_tss_env[c(2,5,6)]) # environmental data: Age, SST, TOC only
otu_tss_t <- as.data.frame(t(otu_tss_env[-(1:6)])) # Transpose abundance matrix taxa x samples
env_t <- as.data.frame(t(otu_tss_env[c(2,5,6)]))
match(rownames(env_scaled), rownames(otu_tss_env))

powers <- c(c(1:20), seq(from = 22, to=30, by=5))

sft <- pickSoftThreshold(otu_tss_env[-(1:6)], powerVector = powers, verbose =5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.60, col="red")  # Typical threshold

net <- blockwiseModules(
  otu_tss_env[-(1:6)],
  power = 6,  # Use the value you chose
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.3,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM",
  verbose = 3
)

mergedColors <- labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)

## Network correlation ----
moduleEigengenes <- net$MEs
moduleTraitCor <- cor(moduleEigengenes, env_scaled, use = "pairwise.complete.obs")
module_order <- paste0("ME", sort(as.integer(gsub("ME", "", rownames(moduleTraitCor)))))
ordered_cor <- moduleTraitCor[module_order, ]
textMatrix <- paste0("r = ", round(moduleTraitCor, 2), "\n",
                     "p = ", signif(moduleTraitPval, 2))


dim(textMatrix) <- dim(moduleTraitCor)
dimnames(textMatrix) <- dimnames(moduleTraitCor)
ordered_text <- textMatrix[module_order, ]
# Correlation of module eigengenes, that represent the dominant pattern of the taxa
# in that module. It is the first PC of the abundance data in that module.
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(otu_tss_env[-(1:6)]))
moduleTraitPval
# Student asymptotic p-value <0.01 only for M0 to lat and temp and M4 to TOC and long


# A : Network Correlation Matrix -----------------------------------------------
dev.off()
pdf("Plots/Fig5a.pdf", width = 15 / 2.54, height = 14 / 2.54)
pheatmap(ordered_cor, 
         cluster_rows = F, 
         cluster_cols = F,
         display_numbers = ordered_text,
         #main = "Module–Trait Correlations",
         legend = F,
         labels_col = c("Age", "Latitude", "Longitude", "UK'37 SST", "TOC"),
         labels_row = paste0(0:4),
         fontsize_col = 8,
         fontsize_number = 8,
         angle_col = 0)
dev.off()

## S1 Centrality measure -------------------------------------------------------

softPower <- 6
all_modules_df <- list()

# Loop over each of the 5 modules
for (mod in unique(net$colors[net$colors != "grey"])) {
  
  module_taxa <- names(net$colors)[net$colors == mod]
  expr_module <- otu_tss_env[, module_taxa, drop = FALSE]
  
  # Adjacency matrix with Spearman correlations
  adj_matrix <- abs(cor(expr_module, method = "spearman"))^softPower
  diag(adj_matrix) <- 0
  
  module_colors <- rep("selectedModule", ncol(expr_module))
  names(module_colors) <- colnames(expr_module)
  
  # kWithin
  k_table <- intramodularConnectivity(adj_matrix, module_colors)
  kWithin <- k_table$kWithin
  names(kWithin) <- rownames(k_table)
  
  # Correlation to TOC
  trait <- env_scaled[rownames(expr_module), "TOC"]
  trait_cor <- cor(expr_module, trait, method = "pearson")
  names(trait_cor) <- colnames(expr_module)
  
  # VIP scores
  vip_scores <- vip(res.sPLS)[, 1]
  vip_scores <- vip_scores[names(vip_scores) %in% module_taxa]
  
  # Combine
  df <- data.frame(
    Taxon = names(kWithin),
    kWithin = kWithin,
    CorTOC = trait_cor[names(kWithin)],
    VIP = vip_scores[names(kWithin)],
    Module = paste0("M", mod)
  )
  
  all_modules_df[[paste0("M", mod)]] <- df
}

plot_df_all <- bind_rows(all_modules_df)

# Plot
ggplot(plot_df_all, aes(x = kWithin, y = CorTOC, size = VIP, fill = Module)) +
  geom_point(alpha = 0.9, shape = 21, color = "black") +
  geom_text_repel(aes(label = Taxon), size = 2.5, max.overlaps = 8, show.legend = FALSE) +
  facet_wrap(~ Module) +
  scale_fill_paletteer_d("MexBrewer::Casita1") +
  theme_bw() +
  labs(
    x = "Intramodular Connectivity (kWithin)",
    y = "Pearson Correlation with TOC",
    size = "VIP Score",
    color = "Module",
    #title = "Subnetwork Connectivity vs. TOC Correlation per Module"
    )
  
ggsave("Plots/S1.pdf", width = 25, height = 20, units = "cm")

### Relative module contributions to each sample ----
head(net$colors) # Taxa that fall in each network
net_colors <- as.data.frame(net$color) %>% rownames_to_column("taxa")

# Aggregate abundance per module
core_lat_order <- otu_tss_env %>%
  dplyr::select(core, latitude) %>%
  distinct() %>%
  mutate(abs_lat = abs(latitude)) %>%
  arrange(abs_lat) %>%
  pull(core)

net_abund <- otu_tss_t %>% 
  rownames_to_column("taxa") %>%
  full_join(net_colors, by = "taxa") %>% 
  rename("module" = "net$color") %>% 
  mutate(module = as.factor(module))

module_contrib <- net_abund %>%
  group_by(module) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  pivot_longer(-module, names_to = "sample", values_to = "contribution") %>% 
  group_by(sample) %>%
  mutate(relative_contribution = contribution / sum(contribution)) %>%  # Normalize in each sample by TSS
  left_join(otu_tss_env[c(1:6)] %>% rownames_to_column("sample"), by = "sample") %>% 
  group_by(core) %>%
  mutate(TOC_scaled = TOC / max(TOC, na.rm = TRUE)) %>% # Normalize in each sample relative to max value
  mutate(Temp_scaled = `UK37_SST_°C` / max(`UK37_SST_°C`, na.rm = TRUE)) %>%  # Normalize in each sample relative to max value
  ungroup() %>% 
  dplyr::mutate(contrib_log = log10((contribution *100)+1)) %>% 
  mutate(core = factor(core, levels = core_lat_order))

ggplot(module_contrib, aes(x = age, y = relative_contribution, fill = module)) +
  #geom_bar(stat = "identity", width = 1) +
  geom_col(position = "stack", color = "black", linewidth = 0.3, width = 1) +
  #geom_line(aes(y = TOC_scaled, color = "TOC"), alpha = 0.3, size = 0.8) +
  #geom_point(aes(y = TOC_scaled, color = "TOC"), alpha = 0.5, size = 1.2, show.legend = F) +
  #geom_line(aes(y = Temp_scaled, color = "UK'37"), alpha = 0.3, size = 0.8) +
  #geom_point(aes(y = Temp_scaled, color = "UK'37"), alpha = 0.5, size = 1.2, show.legend = F) +
  facet_wrap(~ core) +
  labs(#title = "Relative Module Contribution per Sample",
       x = "Sample Age [ka BP]",
       y = "Relative Contribution",
       fill = "Module",
       color = "Parameter") +
  #scale_fill_manual(values = module_colors) +
  paletteer::scale_fill_paletteer_d("MexBrewer::Casita1") +
  #paletteer::scale_fill_paletteer_d("soilpalettes::redox2") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks=seq(2,30,by=4)) +
  #scale_color_manual(values = c("TOC" = "black", "UK'37" = "white")) +
  theme_bw()

ggsave("Plots/Fig5b.pdf", width = 25, height = 15, units = "cm")

# ### B: Module contribution ----
# ggplot(module_contrib, aes(x = age, y = relative_contribution, fill = module)) +
#   geom_bar(stat = "identity") +
#   #geom_line(aes(y = TOC_scaled, color = "TOC"), alpha = 0.3, size = 0.8) +
#   #geom_point(aes(y = TOC_scaled, color = "TOC"), alpha = 0.5, size = 1.2, show.legend = F) +
#   #geom_line(aes(y = Temp_scaled, color = "Temperature"), alpha = 0.3, size = 0.8) +
#   #geom_point(aes(y = Temp_scaled, color = "Temperature"), alpha = 0.5, size = 1.2, show.legend = F) +
#   #facet_wrap(~ core) +
#   facet_grid(module ~ core) +
#   labs(#title = "Relative Module Contribution per Sample", 
#        x = "Sample Age [ky BP]", 
#        y = "Relative Contribution",
#        fill = "Module",
#        color = "Parameter") +
#   #scale_fill_manual(values = module_colors) +
#   paletteer::scale_fill_paletteer_d("MexBrewer::Casita1") +
#   #paletteer::scale_fill_paletteer_d("soilpalettes::redox2") +
#   #scale_color_manual(values = c("TOC" = "black", "Temperature" = "white")) +
#   scale_y_continuous(labels = scales::percent, position = "right") +
#   labs(y = NULL) +
#   theme_bw() +
#   theme(strip.text.y = element_blank(),
#         legend.position = "top")
# 
# ggsave("Plots/Fig5b.pdf", width = 25, height = 15, units = "cm")

## S2 Module contribution vs. Parameters ---------------------------------------
ggplot(module_contrib, aes(x = `UK37_SST_°C`, y = TOC, fill = module,  size = contribution*100)) +
  geom_point(shape = 21, stroke = 0.5, color = "black", alpha = 0.9) +
  paletteer::scale_fill_paletteer_d("MexBrewer::Casita1") +
  facet_wrap(.~module) +
  labs(#title = "Relative Module Contribution per Sample",
    x = "UK'37 SST [°C]",
    y = "TOC [%]",
    fill = "Module",
    size = "Relative\ncontribution\n[%]") +
  theme_bw()

ggsave("Plots/S2.pdf", width = 25, height = 15, units = "cm")

### Investigate trends in network contribution ---- 
## by plotting relative contribution against scaled parameter
ggplot(module_contrib, (aes(x = TOC_scaled, y = relative_contribution, color = module))) + 
  geom_point() +
  geom_smooth(method = "lm", se = F, alpha = 0.8) +
  facet_wrap(~ core) +
  paletteer::scale_color_paletteer_d("MexBrewer::Casita1") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw()

## Some stats to proof that the trends are visually there but not really significant
select <- module_contrib %>% filter(core == "GeoB9064-1", module == 4)
cor.test(select$relative_contribution, select$Temp_scaled, method = "pearson")


# Fig 6: Key Players -----------------------------------------------------------
## Filter taxa based on cor and VIP ----
modules <- setdiff(unique(net$colors), "grey")
plot_df_list <- list()
for (mod in modules) {
  # Step 1: Get taxa in this module
  module_taxa <- names(net$colors)[net$colors == mod]
  expr_module <- otu_tss_env[, module_taxa]
  
  # Step 2: Create adjacency matrix (Spearman + soft-threshold)
  cor_matrix <- cor(expr_module, method = "spearman", use = "pairwise.complete.obs")
  adj_matrix <- abs(cor_matrix)^softPower
  diag(adj_matrix) <- 0
  
  # Step 3: Compute intramodular connectivity (kWithin)
  module_colors <- rep("selectedModule", ncol(expr_module))
  names(module_colors) <- colnames(expr_module)
  k_table <- intramodularConnectivity(adj_matrix, module_colors)
  kWithin <- k_table$kWithin
  names(kWithin) <- rownames(k_table)
  
  # Step 4: Correlation to environmental variable (e.g., TOC)
  trait <- env_scaled[rownames(expr_module), "TOC"]
  
  # Initialize vectors
  trait_cor <- numeric(length = ncol(expr_module))
  p_vals <- numeric(length = ncol(expr_module))
  
  names(trait_cor) <- colnames(expr_module)
  names(p_vals) <- colnames(expr_module)
  
  for (taxon in colnames(expr_module)) {
    res <- cor.test(expr_module[[taxon]], trait, method = "pearson")
    trait_cor[taxon] <- res$estimate
    p_vals[taxon] <- res$p.value
  }
  
  # Step 5: VIP scores for current module taxa (from sPLS)
  vip_scores <- vip(res.sPLS)[, 1]
  vip_scores <- vip_scores[names(vip_scores) %in% module_taxa]
  
  # Step 6: Combine into a data frame
  df_mod <- data.frame(
    Taxon = names(kWithin),
    Module = mod,
    kWithin = kWithin,
    CorTOC = trait_cor[names(kWithin)],
    CorPval = p_vals[names(kWithin)],
    VIP = vip_scores[names(kWithin)]
  )
  
  plot_df_list[[as.character(mod)]] <- df_mod
}
all_taxa <- bind_rows(plot_df_list)

all_taxa <- all_taxa %>%
  mutate(
    is_MAG = grepl("[0-9]", Taxon),
    ecology = ifelse(is_MAG, "unknown", NA)  # assign 'unknown' to MAGs
  )

filtered_taxa <- all_taxa %>%
  filter(
    (CorTOC > 0.2 | CorTOC < -0.2),
    CorPval < 0.05,
    VIP >= 1
  )

# Assign known ecological roles for non-MAG taxa
known_ecology <- data.frame(
  Taxon = c(
    "B_Pelagibacter", "B_Prochlorococcus_A", "B_Prochlorococcus_B", "B_Synechococcus_C", 
    "B_Prochlorococcus_C", "B_Patiriisocius", "B_Geodermatophilus", "B_Kiritimatiellae", 
    "B_Cytophagales", "B_Amylibacter", "B_Mariniblastus", "B_Tateyamaria", 
    "B_Winogradskyella", "B_Alteromonas_E", "B_Cyanobacteria", "B_Planktotalea", 
    "B_Pseudothioglobus", "A_Nitrosopelagicus"
  ),
  ecology = c(
    "heterotroph", "photosynthetic", "photosynthetic", "photosynthetic",
    "photosynthetic", "heterotroph", "heterotroph", "heterotroph",
    "heterotroph", "heterotroph", "heterotroph", "heterotroph",
    "heterotroph", "heterotroph", "photosynthetic", "heterotroph",
    "s_oxidiser", "n_oxidiser"
  )
)

filtered_taxa <- filtered_taxa %>%
  left_join(known_ecology, by = "Taxon") %>% 
  mutate(
    ecology = ifelse(!is.na(ecology.y), ecology.y, ecology.x),
    ecology_known = !is.na(ecology.y)
  ) %>%
  dplyr::select(-ecology.x, -ecology.y)

taxa_presence <- colSums(otu_tss_env[, filtered_taxa$Taxon, drop = FALSE] > 0)

filtered_taxa <- filtered_taxa %>%
  mutate(N_sample = taxa_presence[Taxon])

known_taxa <- filtered_taxa %>% filter(ecology_known == TRUE)
unknown_taxa <- filtered_taxa %>% filter(ecology_known == FALSE)

# Plot unknown taxa first as gray points without labels
ggplot() +
  geom_point(data = unknown_taxa, 
             aes(x = VIP, y = CorTOC, size = N_sample, color = ecology),
             alpha = 0.6) +
  # Plot known taxa colored by ecology with labels
  geom_point(data = known_taxa, 
             aes(x = VIP, y = CorTOC, color = ecology, size = N_sample), alpha = 0.9) +
  geom_text_repel(data = known_taxa, 
                  aes(x = VIP, y = CorTOC, label = Taxon, color = ecology),
                  show.legend = FALSE,
                  size = 3, max.overlaps = 5, box.padding = 0.5) +
  annotate("text", x = 2.455, y = -0.13, label = paste0("Classified Taxa: ", nrow(known_taxa))) +
  annotate("text", x = 2.5, y = -0.2, label = paste0("Unclassified MAG's: ", nrow(unknown_taxa))) +
  #paletteer::scale_color_paletteer_d("rockthemes::husker") +
  scale_color_manual(
    name = "Ecology",
    values = c(
      "heterotroph" = "#422537FF",
      "photosynthetic" = "#86A556FF",
      "s_oxidiser" = "#D58078FF",
      "n_oxidiser" = "#624FB0FF",
      "unknown" = "grey"
    ),
    breaks = c("photosynthetic", "heterotroph", "s_oxidiser", "n_oxidiser", "unknown"),
    labels = c("Photosynthetic", "Heterotroph", "Sulfur Oxidiser", "Nitrite Oxidiser", "Unknown")
  ) +
  theme_bw() +
  labs(
    x = "VIP Score",
    y = "Pearson Correlation to TOC",
    color = "Ecology",
    size = "Occurance\nin samples"
  )

ggsave("Plots/Fig6.pdf", width = 27.75, height = 15, units = "cm")

# Fig 7: Relative Abundance ----------------------------------------------------
otu_tss_env_filtered <- otu_tss_env[, known_taxa$Taxon] %>% 
  rownames_to_column("sample") %>%
  pivot_longer(cols = -sample, names_to = "Taxon", values_to = "rel_abundance") %>% 
  inner_join(otu_tss_env[(1:6)] %>% rownames_to_column("sample")) %>% 
  inner_join(known_ecology) %>% 
  group_by(core) %>%
  mutate(TOC_scaled = TOC / max(TOC, na.rm = TRUE)) %>% # Normalize in each sample relative to max value
  mutate(Temp_scaled = `UK37_SST_°C` / max(`UK37_SST_°C`, na.rm = TRUE)) %>%  # Normalize in each sample relative to max value
  ungroup() %>% 
  mutate(core = fct_reorder(as.factor(core), abs(latitude)))


ggplot(otu_tss_env_filtered, aes(x = age, y = rel_abundance, fill = ecology)) +
  geom_col(position = "stack", color = "black", linewidth = 0.3, width = 1) +
  geom_line(aes(y = TOC_scaled, color = "TOC"), alpha = 0.3, size = 0.8) +
  geom_point(aes(y = TOC_scaled, color = "TOC"), alpha = 0.5, size = 1.2, show.legend = F) +
  geom_line(aes(y = Temp_scaled, color = "UK'37 SST"), alpha = 0.3, size = 0.8) +
  geom_point(aes(y = Temp_scaled, color = "UK'37 SST"), alpha = 0.5, size = 1.2, show.legend = F) +
  theme_bw() +
  xlab("Age [ka BP]") +
  ylab("Proportion") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size=7),
    text = element_text(size = 12)) +
  facet_wrap(~core) +
  #paletteer::scale_fill_paletteer_d("rockthemes::husker") +
  scale_fill_manual(
    name = "Ecology",
    values = c(
      "heterotroph" = "#422537FF",
      "photosynthetic" = "#86A556FF",
      "s_oxidiser" = "#D58078FF",
      "n_oxidiser" = "#624FB0FF",
      "unknown" = "grey"
    ),
    breaks = c("photosynthetic", "heterotroph", "s_oxidiser", "n_oxidiser", "unknown"),
    labels = c("Photosynthetic", "Heterotroph", "Sulfur Oxidiser", "Nitrite Oxidiser", "Unknown")
  ) +
  scale_color_manual(name = "Parameter", values = c("TOC" = "black", "UK'37 SST" = "grey")) +
  #scale_fill_manual(name="ecology", values=paletteTax(length(unique(tax_rel_abund_genus$genus)))) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks=seq(2,30,by=4)) +

ggsave("Plots/Fig7.pdf", width = 25, height = 15, units = "cm")
