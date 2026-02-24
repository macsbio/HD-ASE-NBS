# Requires: dplyr, ggplot2, ggrepel
# install.packages(c("dplyr","ggplot2","ggrepel"))

aseManhattan_top20 <- function(input,
                               statvalue = c("pvalue", "qvalue"),
                               annotation,
                               gws = 0.05,
                               filename = "manhattan",
                               top_n_labels = 20) {
  
  statvalue <- match.arg(statvalue)
  
  # join by rsid (annotation must include rsid, chromosome, position, hgnc_symbol)
  dat <- dplyr::left_join(input, annotation, by = "rsid")
  
  # fallback join by chr+pos if rsid match is poor
  if ("hgnc_symbol" %in% names(dat)) {
    frac_na <- mean(is.na(dat$hgnc_symbol))
    if (isTRUE(frac_na > 0.8) &&
        all(c("chromosome", "position") %in% names(input)) &&
        all(c("chromosome", "position", "hgnc_symbol") %in% names(annotation))) {
      
      dat <- dplyr::left_join(
        input,
        dplyr::distinct(annotation[, c("chromosome", "position", "hgnc_symbol")]),
        by = c("chromosome", "position")
      )
    }
  }
  
  dat <- dat %>%
    dplyr::arrange(chromosome, position) %>%
    dplyr::mutate(
      chromosome = as.character(chromosome),
      chromosome = factor(chromosome, levels = as.character(1:22)),
      position   = as.numeric(position),
      STAT       = as.numeric(.data[[statvalue]]),
      LOGP       = -log10(STAT)
    )
  
  # cumulative bp
  dat$BPcum <- NA_real_
  s <- 0
  for (chr in levels(dat$chromosome)) {
    chr_max <- max(dat$position[dat$chromosome == chr], na.rm = TRUE)
    dat$BPcum[dat$chromosome == chr] <- dat$position[dat$chromosome == chr] + s
    s <- s + chr_max
  }
  
  # top N labels (global)
  lab <- dat %>%
    dplyr::filter(!is.na(hgnc_symbol), is.finite(LOGP), !is.na(STAT)) %>%
    dplyr::arrange(STAT) %>%
    dplyr::slice(1:min(top_n_labels, dplyr::n()))
  
  # plot
  manhattanplot <- ggplot2::ggplot(dat, ggplot2::aes(x = BPcum, y = LOGP, colour = chromosome)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_continuous(breaks = NULL, labels = NULL) +   # remove chromosome labels
    ggplot2::scale_color_manual(values = rep(c("#276EBF", "#183059"), 24)) +
    ggplot2::geom_hline(yintercept = -log10(gws), linetype = 2) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = NULL,
      y = expression(paste("-log"[10], "(p)")),
      title = if ("Analysis" %in% names(dat)) paste0(unique(dat$Analysis)) else NULL
    )
  
  if (nrow(lab) > 0) {
    manhattanplot <- manhattanplot +
      ggrepel::geom_text_repel(
        data = lab,
        ggplot2::aes(label = hgnc_symbol),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.1
      )
  }
  
  grDevices::tiff(paste0(filename, ".tif"), res = 300, width = 175, height = 85, units = "mm")
  print(manhattanplot)
  grDevices::dev.off()
}

# --- Run ---
aseManhattan_top20(
  input = diffAseRes,
  statvalue = "pvalue",
  annotation = annTest,
  gws = 0.05,
  filename = "1_vs_2_Manhattan",
  top_n_labels = 20
)
