# ===============================
# Clinical plots + 3-cluster association tests
# ===============================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# -------------------------------
# 1) Input: metadata table
# -------------------------------
# Assumes you have a data.frame called `metadata`
# with a cluster column called `groups` (values "1","2","3")
# and the clinical variables as numeric columns.


# set working dir to the folder containing this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# paths relative to the script directory
data_dir <- file.path(getwd(), "Data")
out_dir  <- file.path(getwd(), "output")

# create output directory if missing
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# (optional) stop if Data folder is missing
if (!dir.exists(data_dir)) stop("Missing Data folder: ", data_dir)


##import clinical data 
clinical_data <- read.delim(file.path(out_dir, "cleaned-metadata.txt"), sep ="\t")




clinical_data$cluster <- factor(clinical_data$cluster, levels = c("1","2","3"))

# Example variable list (replace with your real column names)
clinical_vars_names <- c(
  "cluster",
  "age_of_death",
  "sex",
  "pmi",
  "RIN",
  "age_of_onset",
  "cag",
  "Duration",
  "h.v_cortical_score",
  "h.v_striatal_score",
  "vonsattel_grade",
  "HDSubject.grade",
  "HDSubject.striatal_score",
  "HDSubject.cortical_score",
  "HDSubject.cag_adj_onset",
  "HDSubject.cag_adj_striatal_score",
  "HDSubject.cag_adj_cortical_score",
  "mRNASeq.RIN"
)

# Keep only variables that exist
clinical_data_sub <- clinical_data[, clinical_vars_names[clinical_vars_names %in% colnames(clinical_data)], drop = FALSE]

clinical_data_sub$cluster <- factor(clinical_data$cluster, levels = c("1","2","3"))
clinical_data_sub$sex <- factor(clinical_data$sex, levels = c("M"))
clinical_data_sub$vonsattel_grade <- factor(clinical_data$vonsattel_grade, levels = c("3", "4"))
clinical_data_sub$HDSubject.grade <- factor(clinical_data$HDSubject.grade, levels = c("3", "4"))

# Numeric variables only (exclude cluster)
num_vars <- names(clinical_data_sub)[sapply(clinical_data_sub, is.numeric)]
num_vars <- setdiff(num_vars, "cluster")

# If you want to drop duplicate versions, uncomment and edit:
num_vars <- setdiff(num_vars, c("HDSubject.striatal_score", "HDSubject.cortical_score", "mRNASeq.RIN"))

clinical_long <- clinical_data_sub %>%
  dplyr::select(cluster, dplyr::all_of(num_vars)) %>%
  tidyr::pivot_longer(
    cols = -cluster,
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::filter(!is.na(value))

# Pretty facet labels (edit as needed)
pretty_labels <- c(
  age_of_death = "Age of Death",
  pmi = "Post-mortem Interval",
  RIN = "RIN",
  mRNASeq.RIN = "RNA Integrity Number",
  age_of_onset = "Age of Onset",
  cag = "CAG repeats",
  Duration = "Disease Duration",
  h.v_cortical_score = "Cortical score",
  h.v_striatal_score = "Striatal score",
  HDSubject.striatal_score = "Striatal score (HDSubject)",
  HDSubject.cortical_score = "Cortical score (HDSubject)",
  HDSubject.cag_adj_onset = "CAG adjusted onset",
  HDSubject.cag_adj_striatal_score = "CAG adjusted striatal score",
  HDSubject.cag_adj_cortical_score = "CAG adjusted cortical score"
)

clinical_long <- clinical_long %>%
  dplyr::mutate(
    variable_label = dplyr::recode(variable, !!!pretty_labels, .default = variable)
  )

# -----------------------------
# 2) Kruskal-Wallis + epsilon²
# -----------------------------
eps2_kw <- function(H, n, k) {
  e <- (as.numeric(H) - k + 1) / (n - k)
  ifelse(is.finite(e) & e > 0, e, 0)
}

min_n_per_group <- 3  # set to 3 if missingness is low

kw_results <- clinical_long %>%
  dplyr::group_by(variable, variable_label) %>%
  dplyr::group_modify(~{
    d <- .x
    tab <- table(d$cluster)
    
    if (length(tab) < 3 || any(tab < min_n_per_group)) {
      return(dplyr::tibble(
        H = NA_real_, p.value = NA_real_, epsilon2 = NA_real_,
        n_total = nrow(d),
        n_c1 = as.integer(ifelse(is.na(tab["1"]), 0, tab["1"])),
        n_c2 = as.integer(ifelse(is.na(tab["2"]), 0, tab["2"])),
        n_c3 = as.integer(ifelse(is.na(tab["3"]), 0, tab["3"]))
      ))
    }
    
    kt <- stats::kruskal.test(value ~ cluster, data = d)
    H <- unname(kt$statistic)
    n <- nrow(d)
    k <- length(tab)
    
    dplyr::tibble(
      H = H,
      p.value = kt$p.value,
      epsilon2 = eps2_kw(H, n, k),
      n_total = n,
      n_c1 = as.integer(unname(tab["1"])),
      n_c2 = as.integer(unname(tab["2"])),
      n_c3 = as.integer(unname(tab["3"]))
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(padj_BH = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(padj_BH, p.value)

# Optional: inspect results
print(kw_results)

# -----------------------------
# 3) Annotation labels for plot (NO FDR shown)
# -----------------------------
p_labels <- clinical_long %>%
  dplyr::group_by(variable, variable_label) %>%
  dplyr::summarise(
    y_pos = max(value, na.rm = TRUE) + 0.10 * diff(range(value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    kw_results %>% dplyr::select(variable, p.value, epsilon2),  # <- no padj_BH
    by = "variable"
  ) %>%
  dplyr::mutate(
    label = ifelse(
      is.na(p.value),
      "KW p = NA",
      paste0("p = ", signif(p.value, 2),
             "\nε² = ", signif(epsilon2, 2))
    )
  )

# -----------------------------
# 4) Plot (3 clusters)
# -----------------------------
p <- ggplot(clinical_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(trim = FALSE, alpha = 0.9, color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.7) +
  geom_text(
    data = p_labels,
    aes(x = 2, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("1" = "#F8766D", "2" = "#00BA38", "3" = "#00BFC4")) +
  labs(x = "Cluster", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

png("clinical_data.png", width = 1500, height = 1000, res = 120)
p
dev.off()

# Example: sex vs cluster
tab_sex <- table(clinical_data_sub$cluster, clinical_data_sub$HDSubject.grade)
fisher.test(tab_sex)

ggplot(clinical_data_sub, aes(x = cluster, fill = vonsattel_grade)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", x = "Cluster", title = "Sex distribution by cluster") +
  theme_bw()