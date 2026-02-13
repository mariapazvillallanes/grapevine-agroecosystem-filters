# ======================================================================
# 04_figures_main_fungi.R
# Purpose:
#   Generate MAIN figures for culturable fungi (paper-ready):
#   - Fig1: Mean class composition (%) by Management × Medium × Organ
#   - Fig2: Exclusive/shared classes (Commercial vs Natural) by Organ × Management
#   - Fig3: Delta richness (Natural - Commercial) by Organ, faceted by Management
#   - Fig4: PCoA (Jaccard) using objects generated in 03_fungi_analysis.R
#
# Inputs:
#   data/processed/Hongos_withSampleID.csv
#   outputs/objects/fungi_meta.rds
#   outputs/objects/fungi_pcoa_jaccard.rds
#
# Outputs (not tracked by git; outputs/ is ignored):
#   outputs/figures/main/Fig5_FungalComposition.pdf
#   outputs/figures/main/Fig6_ExclusiveSharedClasses_Fungi.pdf
#   outputs/figures/main/Fig7_DeltaRichness_Fungi.pdf
#   outputs/figures/main/Fig8_PCoA_Jaccard_Fungi.pdf
#   outputs/tables/Table_S1_DeltaRichness_Fungi.csv
#   outputs/tables/ExclusiveClasses_Detailed_Fungi.csv
# ======================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# ---- Paths ----
infile <- "data/processed/Hongos_withSampleID.csv"

dir.create("outputs/figures/main", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",       recursive = TRUE, showWarnings = FALSE)

# ---- Consistent factor levels ----
organ_levels  <- c("Bark", "Grape", "Leaf", "Soil")
mgmt_levels   <- c("Conventional", "Organic", "Wild")
medium_levels <- c("Commercial", "Natural")

# ---- Read + trim ----
df <- readr::read_delim(infile, delim = ";", show_col_types = FALSE, trim_ws = TRUE) %>%
  mutate(
    Count = as.numeric(Count),
    Organ = str_trim(Organ),
    Management = str_trim(Management),
    Medium = str_trim(Medium),
    FungalClass = str_trim(FungalClass),
    PlantID = str_trim(as.character(PlantID)),
    SampleID = str_trim(as.character(SampleID))
  ) %>%
  mutate(
    Organ = factor(Organ, levels = organ_levels),
    Management = factor(Management, levels = mgmt_levels),
    Medium = factor(Medium, levels = medium_levels)
  )

# ---- Palette used across figures (organs) ----
organ_cols <- c(
  Bark  = "#E69F00",
  Grape = "#56B4E9",
  Leaf  = "#009E73",
  Soil  = "#0072B2"
)

# ======================================================================
# FIGURE 1 — Mean composition (%), averaging plant-level compositions
# ======================================================================

plant_class <- df %>%
  group_by(PlantID, Management, Organ, Medium, FungalClass) %>%
  summarise(Isolates = sum(Count, na.rm = TRUE), .groups = "drop")

plant_totals <- plant_class %>%
  group_by(PlantID, Management, Organ, Medium) %>%
  summarise(Total = sum(Isolates, na.rm = TRUE), .groups = "drop") %>%
  mutate(Present = Total > 0)

plant_comp <- plant_class %>%
  left_join(plant_totals, by = c("PlantID", "Management", "Organ", "Medium")) %>%
  mutate(Freq = ifelse(Present, Isolates / Total, NA_real_))

all_classes <- sort(unique(plant_comp$FungalClass))

plant_comp0 <- plant_totals %>%
  filter(Present) %>%
  select(PlantID, Management, Organ, Medium, Total, Present) %>%
  left_join(
    plant_comp %>% select(PlantID, Management, Organ, Medium, FungalClass, Freq),
    by = c("PlantID", "Management", "Organ", "Medium")
  ) %>%
  group_by(PlantID, Management, Organ, Medium, Total, Present) %>%
  tidyr::complete(FungalClass = all_classes, fill = list(Freq = 0)) %>%
  ungroup()

top_n <- 6
top_classes <- plant_comp0 %>%
  group_by(FungalClass) %>%
  summarise(GlobalMean = mean(Freq), .groups = "drop") %>%
  arrange(desc(GlobalMean)) %>%
  slice_head(n = top_n) %>%
  pull(FungalClass)

plant_comp2 <- plant_comp0 %>%
  mutate(Class_plot = ifelse(FungalClass %in% top_classes, FungalClass, "Other")) %>%
  group_by(PlantID, Management, Organ, Medium, Total, Present, Class_plot) %>%
  summarise(Freq = sum(Freq), .groups = "drop") %>%
  mutate(Class_plot = factor(Class_plot, levels = c(top_classes, "Other")))

plot_df1 <- plant_comp2 %>%
  group_by(Management, Organ, Medium, Class_plot) %>%
  summarise(
    MeanFreq = mean(Freq),
    n_recovered = n_distinct(PlantID),
    .groups = "drop"
  ) %>%
  tidyr::complete(Management, Organ, Medium, Class_plot, fill = list(MeanFreq = 0, n_recovered = 0)) %>%
  mutate(MeanFreq = pmax(0, pmin(MeanFreq, 1)))

panel_counts <- plant_totals %>%
  group_by(Management, Organ, Medium) %>%
  summarise(
    n_plants = n_distinct(PlantID),
    n_recovered = n_distinct(PlantID[Present]),
    .groups = "drop"
  ) %>%
  tidyr::complete(Management, Organ, Medium, fill = list(n_plants = 0, n_recovered = 0)) %>%
  mutate(label = ifelse(n_recovered == 0, "no recovery", paste0("n=", n_recovered, "/", n_plants)))

base_cols <- c("#4E79A7", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1",
               "#F28E2B", "#9C755F", "#BAB0AC", "#86BCB6")
pal_needed <- c(base_cols[seq_len(length(top_classes))], "grey75")
names(pal_needed) <- c(top_classes, "Other")

p_fig1 <- ggplot(plot_df1, aes(x = Organ, y = MeanFreq * 100, fill = Class_plot)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.78) +
  facet_grid(Management ~ Medium) +
  scale_fill_manual(values = pal_needed, drop = TRUE) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Mean composition (%)", fill = "Fungal class") +
  geom_text(
    data = panel_counts,
    aes(x = Organ, y = 98, label = label),
    inherit.aes = FALSE,
    size = 3,
    vjust = 1
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(0.8, "lines"),
    legend.key.height = unit(0.55, "lines")
  )

ggsave("outputs/figures/main/Fig1_FungalComposition.pdf",
       p_fig1, width = 10.2, height = 7.0, units = "in", device = cairo_pdf)

# ======================================================================
# FIGURE 2 — Exclusive/shared classes (Commercial vs Natural)
# ======================================================================

pa <- df %>%
  group_by(Organ, Management, Medium, PlantID, FungalClass) %>%
  summarise(present = sum(Count, na.rm = TRUE) > 0, .groups = "drop") %>%
  group_by(Organ, Management, Medium, FungalClass) %>%
  summarise(present = any(present), .groups = "drop")

wide <- pa %>%
  mutate(group = paste(Organ, Management, sep = " · ")) %>%
  select(group, Organ, Management, FungalClass, Medium, present) %>%
  tidyr::complete(group, Organ, Management, FungalClass, Medium, fill = list(present = FALSE)) %>%
  pivot_wider(names_from = Medium, values_from = present, values_fill = FALSE)

up_counts <- wide %>%
  mutate(category = case_when(
    Commercial & !Natural ~ "Only Commercial",
    Natural & !Commercial ~ "Only Natural",
    Commercial & Natural  ~ "Commercial & Natural",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(category)) %>%
  count(Organ, Management, category, name = "n") %>%
  mutate(
    Organ = factor(Organ, levels = organ_levels),
    Management = factor(Management, levels = mgmt_levels),
    category = factor(category, levels = c("Only Commercial", "Only Natural", "Commercial & Natural"))
  )

p_fig2 <- ggplot(up_counts, aes(x = category, y = n, fill = Organ)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.5) +
  facet_grid(Management ~ Organ, switch = "x") +
  scale_fill_manual(values = organ_cols) +
  scale_x_discrete(labels = c(
    "Only Commercial" = "Commercial",
    "Only Natural" = "Natural",
    "Commercial & Natural" = "Both"
  )) +
  labs(x = NULL, y = "Number of classes") +
  theme_bw(base_size = 12, base_family = "Times New Roman") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(color = "black"),
    strip.placement = "outside",
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

ggsave("outputs/figures/main/Fig2_ExclusiveSharedClasses_Fungi.pdf",
       p_fig2, width = 18, height = 12, units = "cm", device = cairo_pdf)

exclusive_classes_detailed <- wide %>%
  mutate(Category = case_when(
    Commercial & !Natural ~ "Only Commercial",
    Natural & !Commercial ~ "Only Natural",
    Commercial & Natural  ~ "Commercial & Natural",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Category)) %>%
  select(Organ, Management, Category, FungalClass) %>%
  arrange(Management, Organ, Category, FungalClass)

write.csv2(exclusive_classes_detailed, "outputs/tables/ExclusiveClasses_Detailed_Fungi.csv", row.names = FALSE)

# ======================================================================
# FIGURE 3 — Delta richness (Natural - Commercial), paired within PlantID
# ======================================================================

rich_long <- df %>%
  group_by(PlantID, Management, Organ, Medium) %>%
  summarise(
    Richness = n_distinct(FungalClass[Count > 0]),
    .groups = "drop"
  )

rich_df <- rich_long %>%
  pivot_wider(names_from = Medium, values_from = Richness) %>%
  filter(!is.na(Natural) & !is.na(Commercial)) %>%
  mutate(Delta = Natural - Commercial)

stats_df <- rich_df %>%
  group_by(Management, Organ) %>%
  summarise(
    n = n_distinct(PlantID),
    median_delta = median(Delta),
    Q1 = quantile(Delta, 0.25),
    Q3 = quantile(Delta, 0.75),
    IQR = Q3 - Q1,
    all_equal = sd(Delta) == 0,
    p_value = ifelse(all_equal, 1, wilcox.test(Delta, mu = 0, exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

write.table(stats_df,
            file = "outputs/tables/Table_S1_DeltaRichness_Fungi.csv",
            sep = ";", row.names = FALSE, quote = FALSE)

sig_df <- stats_df %>%
  filter(significance != "ns") %>%
  left_join(
    rich_df %>%
      group_by(Management, Organ) %>%
      summarise(y = max(Delta, na.rm = TRUE) + 0.5, .groups = "drop"),
    by = c("Management", "Organ")
  )

p_fig3 <- ggplot(rich_df, aes(Organ, Delta, color = Organ)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "black", fill = NA) +
  geom_jitter(width = 0.12, size = 2.6) +
  facet_wrap(~Management, nrow = 1) +
  scale_color_manual(values = organ_cols) +
  labs(
    y = expression(Delta~"culturable class richness (Natural - Commercial)"),
    x = NULL
  ) +
  geom_text(
    data = sig_df,
    aes(x = Organ, y = y, label = significance),
    inherit.aes = FALSE,
    size = 5,
    family = "Times New Roman"
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave("outputs/figures/main/Fig3_DeltaRichness_Fungi.pdf",
       p_fig3, width = 8.5, height = 4.2, device = cairo_pdf)

# ======================================================================
# FIGURE 4 — PCoA (Jaccard), using saved objects from analysis
# ======================================================================

meta_fun <- readRDS("outputs/objects/fungi_meta.rds")
pcoa_fun <- readRDS("outputs/objects/fungi_pcoa_jaccard.rds")

meta_ord <- meta_fun %>% filter(N > 0)

eig <- pcoa_fun$eig
var1 <- 100 * eig[1] / sum(eig[eig > 0])
var2 <- 100 * eig[2] / sum(eig[eig > 0])

ord <- as.data.frame(pcoa_fun$points)
colnames(ord) <- c("PCoA1", "PCoA2")
ord <- bind_cols(meta_ord, ord)

p_fig4 <- ggplot(ord, aes(PCoA1, PCoA2, color = Organ, shape = Management)) +
  geom_point(size = 3, alpha = 0.55,
             position = position_jitter(width = 0.008, height = 0.008)) +
  scale_color_manual(values = organ_cols) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  labs(
    x = paste0("PCoA1 (", round(var1, 1), "%)"),
    y = paste0("PCoA2 (", round(var2, 1), "%)"),
    color = "Organ",
    shape = "Management"
  )

ggsave("outputs/figures/main/Fig4_PCoA_Jaccard_Fungi.pdf",
       p_fig4, width = 7.2, height = 5.2, device = cairo_pdf)

message("Main fungi figures exported to outputs/figures/main/")

