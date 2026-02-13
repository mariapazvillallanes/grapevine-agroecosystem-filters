# ======================================================================
# 02_bacteria_analysis.R
# Purpose:
#   Perform class-level community analyses for culturable bacteria:
#   - Build sample x class matrix (counts)
#   - Isolation success (N)
#   - Alpha diversity (Richness, Shannon) + LMM (descriptive)
#   - Beta diversity (Jaccard presence/absence; Bray–Curtis on relative abundances)
#   - PERMANOVA with restricted permutations (blocked by PlantUID)
#   - PERMDISP (dispersion checks)
#   - Sensitivity analysis (N >= 5) for Jaccard PERMANOVA
#
# Inputs:
#   data/processed/Bacterias_withSampleID.csv   (semicolon-separated)
#
# Outputs:
#   outputs/objects/bacteria_meta.rds
#   outputs/objects/bacteria_matrix_counts.rds
#   outputs/objects/bacteria_matrix_pa.rds
#   outputs/objects/bacteria_dist_jaccard.rds
#   outputs/objects/bacteria_dist_bray.rds
#   outputs/objects/bacteria_pcoa_jaccard.rds
#   outputs/tables/bacteria_permanova_jaccard.csv
#   outputs/tables/bacteria_permanova_bray.csv
#   outputs/tables/bacteria_permdisp.csv
#   outputs/tables/bacteria_sensitivity_Nge5_jaccard.csv
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
input_file <- "data/processed/Bacterias_withSampleID.csv"

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
required_cols <- c("SampleID", "PlantID", "Management", "Medium", "Organ", "BacterialClass", "Count")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns in input file: ", paste(missing_cols, collapse = ", "))
}

# Type harmonization
df <- df %>%
  mutate(
    SampleID       = as.factor(SampleID),
    PlantID        = as.factor(PlantID),
    Management     = as.factor(Management),
    Medium         = as.factor(Medium),
    Organ          = as.factor(Organ),
    BacterialClass = as.factor(BacterialClass),
    Count          = as.integer(Count)
  )

# Quick summary
message("Unique SampleID: ", n_distinct(df$SampleID))
message("Unique classes:  ", n_distinct(df$BacterialClass))

# ---- 2) PlantUID (unique plant identity within management) ----
df <- df %>%
  mutate(
    PlantUID = as.factor(paste(Management, PlantID, sep = "_"))
  )

# ---- 3) Build sample x class count matrix ----
mat_counts <- df %>%
  select(SampleID, BacterialClass, Count) %>%
  group_by(SampleID, BacterialClass) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = BacterialClass, values_from = Count, values_fill = 0) %>%
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
saveRDS(meta,      "outputs/objects/bacteria_meta.rds")
saveRDS(X_counts,  "outputs/objects/bacteria_matrix_counts.rds")

# ---- 4) Alpha diversity (descriptive; strongly driven by N) ----
meta$Richness <- vegan::specnumber(X_counts)
meta$Shannon  <- vegan::diversity(X_counts, index = "shannon")

# LMM checks (descriptive)
m_rich <- lmer(Richness ~ Organ + Management + Medium + scale(N) + (1 | PlantUID), data = meta)
m_shan <- lmer(Shannon  ~ Organ + Management + Medium + scale(N) + (1 | PlantUID), data = meta)

tab_rich <- broom.mixed::tidy(m_rich, effects = "fixed")
tab_shan <- broom.mixed::tidy(m_shan, effects = "fixed")

readr::write_delim(tab_rich, "outputs/tables/bacteria_LMM_richness.csv", delim = ";")
readr::write_delim(tab_shan, "outputs/tables/bacteria_LMM_shannon.csv",  delim = ";")

# ---- 5) Beta diversity ----

# 5.1 Jaccard (presence/absence); exclude N=0 samples
keep_nonzero <- meta$N > 0
X_pa <- (X_counts > 0)[keep_nonzero, , drop = FALSE]
meta_pa <- meta[keep_nonzero, , drop = FALSE]

dist_jacc <- vegan::vegdist(X_pa, method = "jaccard")
saveRDS(X_pa,      "outputs/objects/bacteria_matrix_pa.rds")
saveRDS(dist_jacc, "outputs/objects/bacteria_dist_jaccard.rds")

# PCoA (save object; figure generated elsewhere)
pcoa_jacc <- cmdscale(dist_jacc, k = 2, eig = TRUE)
saveRDS(pcoa_jacc, "outputs/objects/bacteria_pcoa_jaccard.rds")

# 5.2 Bray–Curtis on relative abundances (counts -> proportions per sample)
X_rel <- sweep(X_counts, 1, rowSums(X_counts), FUN = "/")
X_rel[is.na(X_rel)] <- 0  # handle N=0 samples

# Exclude N=0 for distances
X_rel_nz <- X_rel[keep_nonzero, , drop = FALSE]
meta_rel <- meta[keep_nonzero, , drop = FALSE]

dist_bray <- vegan::vegdist(X_rel_nz, method = "bray")
saveRDS(dist_bray, "outputs/objects/bacteria_dist_bray.rds")

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

# Export PERMANOVA tables
readr::write_delim(
  as.data.frame(adonis_jacc),
  "outputs/tables/bacteria_permanova_jaccard.csv",
  delim = ";"
)

readr::write_delim(
  as.data.frame(adonis_bray),
  "outputs/tables/bacteria_permanova_bray.csv",
  delim = ";"
)

# ---- 7) PERMDISP (dispersion) ----
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

readr::write_delim(permdisp_tbl, "outputs/tables/bacteria_permdisp.csv", delim = ";")

# ---- 8) Sensitivity analysis: N >= 5 (Jaccard) ----
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

readr::write_delim(
  as.data.frame(sens),
  "outputs/tables/bacteria_sensitivity_Nge5_jaccard.csv",
  delim = ";"
)

message("Bacteria analysis completed successfully.")
