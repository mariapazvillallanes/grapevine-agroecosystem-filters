# ======================================================================
# 03_fungi_analysis.R
# Purpose:
#   Perform class-level community analyses for culturable fungi:
#   - Build sample x class matrix (counts)
#   - Isolation success (N)
#   - Alpha diversity (Richness, Shannon) + LMM (descriptive)
#   - Beta diversity (Jaccard presence/absence; Bray–Curtis on relative abundances)
#   - PERMANOVA with restricted permutations (blocked by PlantUID)
#   - PERMDISP (dispersion checks)
#   - Sensitivity analysis (N >= 5) for Jaccard PERMANOVA
#
# Inputs:
#   data/processed/Hongos_withSampleID.csv   (semicolon-separated)
#
# Outputs:
#   outputs/objects/fungi_meta.rds
#   outputs/objects/fungi_matrix_counts.rds
#   outputs/objects/fungi_matrix_pa.rds
#   outputs/objects/fungi_dist_jaccard.rds
#   outputs/objects/fungi_dist_bray.rds
#   outputs/objects/fungi_pcoa_jaccard.rds
#   outputs/tables/fungi_LMM_richness.csv
#   outputs/tables/fungi_LMM_shannon.csv
#   outputs/tables/fungi_permanova_jaccard.csv
#   outputs/tables/fungi_permanova_bray.csv
#   outputs/tables/fungi_permdisp.csv
#   outputs/tables/fungi_sensitivity_Nge5_jaccard.csv
# ======================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(lme4)
  library(lmerTest)
  library(permute)
  library(broom.mixed)
})

# ---- Paths ----
input_file <- "data/processed/Hongos_withSampleID.csv"

dir.create("outputs/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)

# ---- 1) Load data ----
message("Reading: ", input_file)

df <- readr::read_delim(
  file = input_file,
  delim = ";",
  show_col_types = FALSE,
  trim_ws = TRUE
)

# Required columns check
required_cols <- c("SampleID", "PlantID", "Management", "Medium", "Organ", "FungalClass", "Count")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns in input file: ", paste(missing_cols, collapse = ", "))
}

# Type harmonization
df <- df %>%
  mutate(
    SampleID    = as.factor(SampleID),
    PlantID     = as.factor(PlantID),
    Management  = as.factor(Management),
    Medium      = as.factor(Medium),
    Organ       = as.factor(Organ),
    FungalClass = as.factor(FungalClass),
    Count       = as.integer(Count)
  )

message("Unique SampleID: ", n_distinct(df$SampleID))
message("Unique classes:  ", n_distinct(df$FungalClass))

# ---- 2) PlantUID (unique plant identity within management) ----
df <- df %>%
  mutate(
    PlantUID = as.factor(paste(Management, PlantID, sep = "_"))
  )

# ---- 3) Build sample x class count matrix ----
mat_counts <- df %>%
  select(SampleID, FungalClass, Count) %>%
  group_by(SampleID, FungalClass) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = FungalClass, values_from = Count, values_fill = 0) %>%
  arrange(SampleID)

meta <- df %>%
  select(SampleID, PlantUID, PlantID, Management, Medium, Organ) %>%
  distinct() %>%
  arrange(SampleID)

stopifnot(all(meta$SampleID == mat_counts$SampleID))

X_counts <- mat_counts %>%
  select(-SampleID) %>%
  as.data.frame() %>%
  as.matrix()

# Isolation success (total isolates per sample)
meta$N <- rowSums(X_counts)

# Save core objects
saveRDS(meta,     "outputs/objects/fungi_meta.rds")
saveRDS(X_counts, "outputs/objects/fungi_matrix_counts.rds")

# ---- 4) Alpha diversity (descriptive; constrained by N) ----
meta$Richness <- vegan::specnumber(X_counts)
meta$Shannon  <- vegan::diversity(X_counts, index = "shannon")

m_rich <- lmer(Richness ~ Organ + Management + Medium + scale(N) + (1 | PlantUID), data = meta)
m_shan <- lmer(Shannon  ~ Organ + Management + Medium + scale(N) + (1 | PlantUID), data = meta)

tab_rich <- broom.mixed::tidy(m_rich, effects = "fixed")
tab_shan <- broom.mixed::tidy(m_shan, effects = "fixed")

readr::write_delim(tab_rich, "outputs/tables/fungi_LMM_richness.csv", delim = ";")
readr::write_delim(tab_shan, "outputs/tables/fungi_LMM_shannon.csv",  delim = ";")

# ---- 5) Beta diversity ----
keep_nonzero <- meta$N > 0

# 5.1 Jaccard (presence/absence)
X_pa   <- (X_counts > 0)[keep_nonzero, , drop = FALSE]
meta_pa <- meta[keep_nonzero, , drop = FALSE]

dist_jacc <- vegan::vegdist(X_pa, method = "jaccard")

saveRDS(X_pa,      "outputs/objects/fungi_matrix_pa.rds")
saveRDS(dist_jacc, "outputs/objects/fungi_dist_jaccard.rds")

# PCoA object (figure generated in figures script)
pcoa_jacc <- cmdscale(dist_jacc, k = 2, eig = TRUE)
saveRDS(pcoa_jacc, "outputs/objects/fungi_pcoa_jaccard.rds")

# 5.2 Bray–Curtis on relative abundances
X_rel <- sweep(X_counts, 1, rowSums(X_counts), FUN = "/")
X_rel[is.na(X_rel)] <- 0

X_rel_nz <- X_rel[keep_nonzero, , drop = FALSE]
meta_rel <- meta[keep_nonzero, , drop = FALSE]

dist_bray <- vegan::vegdist(X_rel_nz, method = "bray")
saveRDS(dist_bray, "outputs/objects/fungi_dist_bray.rds")

# ---- 6) PERMANOVA with restricted permutations (blocked by PlantUID) ----
perm_j <- how(nperm = 999)
setBlocks(perm_j) <- meta_pa$PlantUID

adonis_jacc <- adonis2(
  dist_jacc ~ Management * Organ + Medium,
  data = meta_pa,
  permutations = perm_j,
  by = "margin"
)

perm_b <- how(nperm = 999)
setBlocks(perm_b) <- meta_rel$PlantUID

adonis_bray <- adonis2(
  dist_bray ~ Management * Organ + Medium,
  data = meta_rel,
  permutations = perm_b,
  by = "margin"
)

readr::write_delim(as.data.frame(adonis_jacc), "outputs/tables/fungi_permanova_jaccard.csv", delim = ";")
readr::write_delim(as.data.frame(adonis_bray), "outputs/tables/fungi_permanova_bray.csv",    delim = ";")

# ---- 7) PERMDISP ----
set.seed(123)

bd_med_j <- betadisper(dist_jacc, meta_pa$Medium)
bd_int_j <- betadisper(dist_jacc, interaction(meta_pa$Management, meta_pa$Organ))
pt_med_j <- permutest(bd_med_j, permutations = 999)
pt_int_j <- permutest(bd_int_j, permutations = 999)

bd_med_b <- betadisper(dist_bray, meta_rel$Medium)
bd_int_b <- betadisper(dist_bray, interaction(meta_rel$Management, meta_rel$Organ))
pt_med_b <- permutest(bd_med_b, permutations = 999)
pt_int_b <- permutest(bd_int_b, permutations = 999)

permdisp_tbl <- tibble(
  Distance = c("Jaccard", "Jaccard", "Bray-Curtis", "Bray-Curtis"),
  Factor   = c("Medium", "Management x Organ", "Medium", "Management x Organ"),
  F        = c(pt_med_j$tab[1,"F"], pt_int_j$tab[1,"F"], pt_med_b$tab[1,"F"], pt_int_b$tab[1,"F"]) |> as.numeric(),
  p_value  = c(pt_med_j$tab[1,"Pr(>F)"], pt_int_j$tab[1,"Pr(>F)"], pt_med_b$tab[1,"Pr(>F)"], pt_int_b$tab[1,"Pr(>F)"]) |> as.numeric()
) %>%
  mutate(
    F = round(F, 3),
    p_value = signif(p_value, 3)
  )

readr::write_delim(permdisp_tbl, "outputs/tables/fungi_permdisp.csv", delim = ";")

# ---- 8) Sensitivity: N >= 5 (Jaccard) ----
keep5 <- meta$N >= 5
X_pa_5 <- (X_counts > 0)[keep5, , drop = FALSE]
meta_5 <- meta[keep5, , drop = FALSE]

dist_jacc_5 <- vegan::vegdist(X_pa_5, method = "jaccard")

perm5 <- how(nperm = 999)
setBlocks(perm5) <- meta_5$PlantUID

sens <- adonis2(
  dist_jacc_5 ~ Management * Organ + Medium,
  data = meta_5,
  permutations = perm5,
  by = "margin"
)

readr::write_delim(as.data.frame(sens), "outputs/tables/fungi_sensitivity_Nge5_jaccard.csv", delim = ";")

message("Fungi analysis completed successfully.")
