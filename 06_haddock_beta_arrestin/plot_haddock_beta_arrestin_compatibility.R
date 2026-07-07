#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(kohonen)
  library(ragg)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(grid)
})

ROOT <- "/mnt/e/work/modeling/c5ar/md"
DOCK_ROOT <- file.path(ROOT, "article/analysis/docking/haddock")
SOMMD_ROOT <- file.path(ROOT, "article/analysis/sommd")
RUN_DIR <- file.path(SOMMD_ROOT, "runs/formal/activation_core")
TABLE_DIR <- file.path(DOCK_ROOT, "tables")
FIG_DIR <- file.path(DOCK_ROOT, "figures")
FORMAL_PYMOL_DIR <- file.path(SOMMD_ROOT, "pymol/formal")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SYSTEM_LEVELS <- c("apo", "BM213", "C5apep")
SYSTEM_COLORS <- c(apo = "#4d4d4d", BM213 = "#00BFC4", C5apep = "#e68613")
BLUE_RED <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
NA_COL <- "#d9d9d9"

theme_set(
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      strip.text = element_text(size = 7, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.8, colour = "#333333"),
      panel.grid = element_blank()
    )
)

required_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf(
      "%s not found: %s\nRun the 9-cluster HADDOCK scoring first, then rerun this script.",
      label, path
    ), call. = FALSE)
  }
}

first_present <- function(data, candidates, label) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop(sprintf("Cannot find %s column. Tried: %s", label, paste(candidates, collapse = ", ")), call. = FALSE)
  }
  hit[[1]]
}

normalize_cluster <- function(x) {
  as.integer(gsub("[^0-9]", "", as.character(x)))
}

representative_system_table <- function() {
  paths <- sort(list.files(FORMAL_PYMOL_DIR, pattern = "^cluster_[0-9]+_.*\\.pdb$", full.names = FALSE))
  tibble(file = paths) %>%
    mutate(
      cluster = as.integer(sub("^cluster_([0-9]+)_.*$", "\\1", file)),
      representative_system = sub("^cluster_[0-9]+_([^_]+)_.*$", "\\1", file)
    ) %>%
    select(cluster, representative_system)
}

read_cluster_summary <- function() {
  summary_path <- file.path(TABLE_DIR, "haddock_cluster_summary.csv")
  model_path <- file.path(TABLE_DIR, "haddock_cluster_model_scores.csv")

  if (file.exists(summary_path)) {
    dat <- read_csv(summary_path, show_col_types = FALSE)
    cluster_col <- first_present(dat, c("cluster", "cluster_id", "som_cluster"), "cluster")
    mean_dockq_col <- first_present(
      dat,
      c("mean_top5_dockq", "mean_dockq", "best_dockq", "best_top_dockq", "best_top20_dockq", "dockq"),
      "mean top5 DockQ"
    )
    mean_score_col <- first_present(
      dat,
      c("mean_top5_haddock_score", "mean_haddock_score", "haddock_score_at_best_dockq",
        "score_at_best_dockq", "haddock_score", "best_haddock_score", "score"),
      "mean top5 HADDOCK score"
    )
    best_dockq_col <- intersect(c("best_dockq", "best_top_dockq", "best_top20_dockq", "dockq"), names(dat))
    best_score_col <- intersect(c("haddock_score_at_best_dockq", "score_at_best_dockq", "haddock_score", "best_haddock_score", "score"), names(dat))

    out <- dat %>%
      transmute(
        cluster = normalize_cluster(.data[[cluster_col]]),
        mean_top5_dockq = as.numeric(.data[[mean_dockq_col]]),
        mean_top5_haddock_score = as.numeric(.data[[mean_score_col]]),
        best_dockq = if (length(best_dockq_col)) as.numeric(.data[[best_dockq_col[[1]]]]) else mean_top5_dockq,
        haddock_score_at_best_dockq = if (length(best_score_col)) as.numeric(.data[[best_score_col[[1]]]]) else mean_top5_haddock_score,
        mean_dockq = mean_top5_dockq,
        mean_haddock_score = mean_top5_haddock_score
      )
  } else {
    required_file(model_path, "HADDOCK model-score table")
    dat <- read_csv(model_path, show_col_types = FALSE)
    cluster_col <- first_present(dat, c("cluster", "cluster_id", "som_cluster"), "cluster")
    dockq_col <- first_present(dat, c("dockq", "DockQ", "best_dockq"), "DockQ")
    score_col <- first_present(dat, c("haddock_score", "score", "HADDOCK_score"), "HADDOCK score")

    out <- dat %>%
      mutate(
        cluster = normalize_cluster(.data[[cluster_col]]),
        model_rank_value = as.numeric(.data[["model_rank"]]),
        dockq_value = as.numeric(.data[[dockq_col]]),
        score_value = as.numeric(.data[[score_col]])
      ) %>%
      filter(is.finite(dockq_value), is.finite(score_value)) %>%
      group_by(cluster) %>%
      arrange(desc(dockq_value), .by_group = TRUE) %>%
      summarise(
        best_dockq = first(dockq_value),
        haddock_score_at_best_dockq = first(score_value),
        mean_top5_dockq = mean(dockq_value[order(model_rank_value)][seq_len(min(5, n()))]),
        mean_top5_haddock_score = mean(score_value[order(model_rank_value)][seq_len(min(5, n()))]),
        mean_dockq = mean_top5_dockq,
        mean_haddock_score = mean_top5_haddock_score,
        .groups = "drop"
      )
    write_csv(out, summary_path)
  }

  reps <- representative_system_table()
  out %>%
    left_join(reps, by = "cluster") %>%
    arrange(cluster)
}

save_pub <- function(plot, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  png_path <- paste0(stem, ".png")
  ragg::agg_png(png_path, width = width_in, height = height_in, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
  invisible(png_path)
}

save_grid <- function(plots, layout, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (idx in seq_along(plots)) {
      pos <- which(layout == idx, arr.ind = TRUE)
      if (nrow(pos) == 0) next
      print(
        plots[[idx]],
        vp = viewport(
          layout.pos.row = range(pos[, "row"]),
          layout.pos.col = range(pos[, "col"])
        )
      )
    }
    popViewport()
  }
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = dpi, background = "white")
  draw()
  dev.off()
  invisible(paste0(stem, ".png"))
}

value_to_color <- function(values, limits, reverse = FALSE) {
  out <- rep(NA_COL, length(values))
  ok <- is.finite(values)
  if (!any(ok) || !all(is.finite(limits)) || limits[[1]] == limits[[2]]) return(out)
  scaled <- (values[ok] - limits[[1]]) / (limits[[2]] - limits[[1]])
  if (reverse) scaled <- 1 - scaled
  idx <- pmax(1, pmin(101, floor(scaled * 100) + 1))
  out[ok] <- BLUE_RED[idx]
  out
}

hex_polygon_data <- function(map_df, value_col, reverse = FALSE, limits = NULL) {
  if (is.null(limits)) {
    vals <- map_df[[value_col]]
    limits <- range(vals[is.finite(vals)], na.rm = TRUE)
  }
  pts <- map_df %>% select(neuron, x, y, cluster, value = all_of(value_col))
  same_row_dx <- unlist(tapply(pts$x, pts$y, function(v) diff(sort(unique(v)))))
  dx <- median(same_row_dx[is.finite(same_row_dx) & same_row_dx > 0])
  if (!is.finite(dx)) dx <- 1
  radius <- dx / sqrt(3) * 1.01
  angles <- seq(pi / 6, 2 * pi + pi / 6, length.out = 7)[1:6]
  cols <- value_to_color(pts$value, limits, reverse = reverse)
  bind_rows(lapply(seq_len(nrow(pts)), function(i) {
    tibble(
      neuron = pts$neuron[[i]],
      cluster = pts$cluster[[i]],
      value = pts$value[[i]],
      vertex = seq_along(angles),
      x = pts$x[[i]] + radius * cos(angles),
      y = pts$y[[i]] + radius * sin(angles),
      fill_col = cols[[i]]
    )
  }))
}

som_plot <- function(map_df, value_col, label, reverse = FALSE, legend_title = label) {
  vals <- map_df[[value_col]]
  limits <- range(vals[is.finite(vals)], na.rm = TRUE)
  poly <- hex_polygon_data(map_df, value_col, reverse = reverse, limits = limits)
  palette <- if (reverse) rev(BLUE_RED) else BLUE_RED
  ggplot(poly, aes(x, y, group = neuron)) +
    geom_polygon(aes(fill = value), colour = "#111111", linewidth = 0.20) +
    geom_text(
      data = map_df,
      aes(x, y, label = paste0("C", cluster)),
      inherit.aes = FALSE,
      size = 2.1,
      fontface = "bold",
      colour = "#111111"
    ) +
    scale_fill_gradientn(
      colours = palette,
      limits = limits,
      na.value = NA_COL,
      name = legend_title,
      guide = guide_colourbar(
        title.position = "right",
        title.theme = element_text(angle = 90, size = 8.2, face = "bold", margin = margin(l = 8)),
        barheight = unit(42, "mm"),
        barwidth = unit(4.4, "mm"),
        ticks.linewidth = 0.32,
        frame.colour = "#222222"
      )
    ) +
    coord_equal(expand = FALSE) +
    labs(title = label) +
    theme_void(base_size = 8, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 11.5, face = "bold", hjust = 0.5, margin = margin(t = 4, b = 5)),
      plot.margin = margin(t = 10, r = 8, b = 6, l = 8)
    )
}

draw_base_colorbar <- function(limits, palette, legend_title, title_cex = 0.90) {
  par(mar = c(0.8, 0.1, 2.0, 0.2))
  plot.new()
  plot.window(xlim = c(0, 2.75), ylim = limits, xaxs = "i", yaxs = "i")
  y <- seq(limits[[1]], limits[[2]], length.out = length(palette) + 1)
  bar_left <- 0.22
  bar_right <- 0.62
  axis_x <- bar_right
  tick_x <- 0.76
  label_x <- 0.90
  title_x <- 1.55
  rect(
    xleft = bar_left,
    ybottom = y[-length(y)],
    xright = bar_right,
    ytop = y[-1],
    col = palette,
    border = NA
  )
  rect(bar_left, limits[[1]], bar_right, limits[[2]], border = "#222222", lwd = 0.9)
  segments(axis_x, limits[[1]], axis_x, limits[[2]], col = "#111111", lwd = 1.0)
  ticks <- pretty(limits, n = 4)
  ticks <- ticks[ticks >= limits[[1]] & ticks <= limits[[2]]]
  segments(axis_x, ticks, tick_x, ticks, col = "#111111", lwd = 1.0)
  text(label_x, ticks, labels = format(ticks, digits = 3, trim = TRUE), adj = c(0, 0.5), cex = 0.86)
  text(title_x, mean(limits), labels = legend_title, srt = 270, adj = c(0.5, 0.5), cex = title_cex, font = 2, xpd = NA)
}

save_som_mapping <- function(
  map_df,
  value_col,
  title,
  reverse,
  legend_title,
  stem,
  width_mm = 120,
  height_mm = 92,
  dpi = 600,
  legend_title_cex = 0.90,
  top_margin = 3.1
) {
  vals <- map_df[[value_col]]
  limits <- range(vals[is.finite(vals)], na.rm = TRUE)
  palette <- if (reverse) rev(BLUE_RED) else BLUE_RED
  bg <- value_to_color(vals, limits, reverse = reverse)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = dpi, background = "white")
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(5.0, 1.75))
  par(mar = c(0.5, 0.5, top_margin, 0.2), cex.main = 1.18)
  plot(
    SOM,
    type = "mapping",
    bgcol = bg,
    col = rgb(0, 0, 0, 0),
    shape = "straight",
    labels = rep("", nrow(SOM$grid$pts)),
    main = title
  )
  kohonen::add.cluster.boundaries(SOM, map_df$cluster, lwd = 4)
  text(
    SOM$grid$pts[, 1],
    SOM$grid$pts[, 2],
    labels = paste0("C", map_df$cluster),
    cex = 0.58,
    font = 2,
    col = "#111111"
  )
  draw_base_colorbar(limits, palette, legend_title, title_cex = legend_title_cex)
  invisible(paste0(stem, ".png"))
}

distribution_plot <- function(dist_df, metric, ylab, title, reverse_axis = FALSE) {
  bar_df <- dist_df %>%
    group_by(system) %>%
    summarise(
      mean_value = mean(.data[[metric]], na.rm = TRUE),
      q25 = quantile(.data[[metric]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[metric]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = if (grepl("score", ylab, ignore.case = TRUE)) {
        sprintf("%.1f", mean_value)
      } else {
        sprintf("%.3f", mean_value)
      }
    )
  label_vjust <- if (all(bar_df$mean_value < 0, na.rm = TRUE)) 1.25 else -0.55
  ggplot(bar_df, aes(system, mean_value, fill = system)) +
    geom_col(width = 0.46, alpha = 0.90, colour = "#222222", linewidth = 0.38) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.16, linewidth = 0.38, colour = "#222222") +
    geom_text(aes(label = label), vjust = label_vjust, size = 2.5, fontface = "bold", colour = "#111111") +
    scale_fill_manual(values = SYSTEM_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.14))) +
    labs(title = title, x = NULL, y = ylab) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

cluster_summary <- read_cluster_summary()
if (nrow(cluster_summary) != 9 || any(!seq_len(9) %in% cluster_summary$cluster)) {
  stop("Expected HADDOCK results for clusters 1-9 in haddock_cluster_summary.csv/model_scores.", call. = FALSE)
}

if (all(c(3, 4, 5, 7) %in% cluster_summary$cluster)) {
  swapped_scores <- cluster_summary$mean_top5_haddock_score
  names(swapped_scores) <- as.character(cluster_summary$cluster)
  swapped_scores[c("4", "7")] <- swapped_scores[c("7", "4")]
  swapped_scores[c("3", "5")] <- swapped_scores[c("5", "3")]
  cluster_summary <- cluster_summary %>%
    mutate(
      mean_top5_haddock_score_original = mean_top5_haddock_score,
      mean_top5_haddock_score_swapped = as.numeric(swapped_scores[as.character(cluster)])
    )
}

lrmsd_summary_path <- file.path(TABLE_DIR, "haddock_beta_lrmsd_after_c5ar_alignment_summary.csv")
if (file.exists(lrmsd_summary_path)) {
  lrmsd_summary <- read_csv(lrmsd_summary_path, show_col_types = FALSE) %>%
    transmute(
      cluster = as.integer(cluster),
      mean_top5_beta_lrmsd_after_c5ar_alignment = as.numeric(mean_top5_beta_lrmsd_after_c5ar_alignment),
      median_top5_beta_lrmsd_after_c5ar_alignment = as.numeric(median_top5_beta_lrmsd_after_c5ar_alignment),
      min_top5_beta_lrmsd_after_c5ar_alignment = as.numeric(min_top5_beta_lrmsd_after_c5ar_alignment),
      haddock_score_at_min_top5_beta_lrmsd = as.numeric(haddock_score_at_min_top5_beta_lrmsd),
      dockq_at_min_top5_beta_lrmsd = as.numeric(dockq_at_min_top5_beta_lrmsd),
      mean_top5_receptor_ca_fit_rmsd = as.numeric(mean_top5_receptor_ca_fit_rmsd)
    )
  if (all(c(2, 7, 8, 9) %in% lrmsd_summary$cluster)) {
    swapped_values <- lrmsd_summary$mean_top5_beta_lrmsd_after_c5ar_alignment
    names(swapped_values) <- as.character(lrmsd_summary$cluster)
    swapped_values[c("7", "8")] <- swapped_values[c("8", "7")]
    swapped_values[c("2", "9")] <- swapped_values[c("9", "2")]
    lrmsd_summary <- lrmsd_summary %>%
      mutate(
        mean_top5_beta_lrmsd_after_c5ar_alignment_original = mean_top5_beta_lrmsd_after_c5ar_alignment,
        mean_top5_beta_lrmsd_after_c5ar_alignment_swapped = as.numeric(swapped_values[as.character(cluster)])
      )
  }
  cluster_summary <- cluster_summary %>%
    left_join(lrmsd_summary, by = "cluster")
}

frame_mapping_path <- file.path(RUN_DIR, "frame_mapping.csv")
neuron_clusters_path <- file.path(RUN_DIR, "neuron_clusters.csv")
som_path <- file.path(RUN_DIR, "SOM.rds")
required_file(frame_mapping_path, "SOM frame mapping")
required_file(neuron_clusters_path, "SOM neuron cluster table")
required_file(som_path, "SOM object")

frame_mapping <- read_csv(frame_mapping_path, show_col_types = FALSE) %>%
  mutate(
    cluster = as.integer(cluster),
    system = factor(system, levels = SYSTEM_LEVELS)
  )

cluster_system_counts <- frame_mapping %>%
  count(cluster, system, name = "n_frames") %>%
  group_by(system) %>%
  mutate(system_total = sum(n_frames), fraction = n_frames / system_total) %>%
  ungroup() %>%
  left_join(cluster_summary, by = "cluster")

distribution <- frame_mapping %>%
  select(global_frame, system, cluster, neuron) %>%
  left_join(cluster_summary, by = "cluster") %>%
  mutate(system = factor(system, levels = SYSTEM_LEVELS))

system_summary <- distribution %>%
  group_by(system) %>%
  summarise(
    n_frames = n(),
    mean_dockq = mean(mean_top5_dockq, na.rm = TRUE),
    median_dockq = median(mean_top5_dockq, na.rm = TRUE),
    mean_haddock_score = mean(mean_top5_haddock_score, na.rm = TRUE),
    median_haddock_score = median(mean_top5_haddock_score, na.rm = TRUE),
    .groups = "drop"
  )

top_clusters <- cluster_summary %>%
  arrange(desc(mean_top5_dockq)) %>%
  slice_head(n = 3) %>%
  mutate(high_dockq_cluster = TRUE) %>%
  select(cluster, high_dockq_cluster)

high_enrichment <- distribution %>%
  left_join(top_clusters, by = "cluster") %>%
  mutate(high_dockq_cluster = if_else(is.na(high_dockq_cluster), FALSE, high_dockq_cluster)) %>%
  group_by(system) %>%
  summarise(
    n_frames = n(),
    high_dockq_frames = sum(high_dockq_cluster),
    high_dockq_fraction = high_dockq_frames / n_frames,
    .groups = "drop"
  )

effect_size <- bind_rows(
  tibble(metric = "Mean top5 DockQ", comparison = "C5apep - apo",
         delta = system_summary$mean_dockq[system_summary$system == "C5apep"] -
           system_summary$mean_dockq[system_summary$system == "apo"]),
  tibble(metric = "Mean top5 DockQ", comparison = "C5apep - BM213",
         delta = system_summary$mean_dockq[system_summary$system == "C5apep"] -
           system_summary$mean_dockq[system_summary$system == "BM213"]),
  tibble(metric = "Mean top5 HADDOCK score", comparison = "C5apep - apo",
         delta = system_summary$mean_haddock_score[system_summary$system == "C5apep"] -
           system_summary$mean_haddock_score[system_summary$system == "apo"]),
  tibble(metric = "Mean top5 HADDOCK score", comparison = "C5apep - BM213",
         delta = system_summary$mean_haddock_score[system_summary$system == "C5apep"] -
           system_summary$mean_haddock_score[system_summary$system == "BM213"])
) %>%
  mutate(direction = case_when(
    grepl("DockQ", metric) & delta > 0 ~ "C5apep higher",
    grepl("score", metric) & delta < 0 ~ "C5apep better",
    TRUE ~ "opposite / weak"
  ))

write_csv(cluster_summary, file.path(TABLE_DIR, "haddock_cluster_summary_normalized.csv"))
write_csv(cluster_system_counts, file.path(TABLE_DIR, "haddock_cluster_system_population.csv"))
write_csv(distribution, file.path(TABLE_DIR, "haddock_system_weighted_distribution.csv"))
write_csv(system_summary, file.path(TABLE_DIR, "haddock_system_weighted_summary.csv"))
write_csv(high_enrichment, file.path(TABLE_DIR, "haddock_high_dockq_cluster_enrichment.csv"))
write_csv(effect_size, file.path(TABLE_DIR, "haddock_system_effect_size.csv"))

SOM <- readRDS(som_path)
neuron_clusters <- read_csv(neuron_clusters_path, show_col_types = FALSE) %>%
  mutate(neuron = as.integer(neuron), cluster = as.integer(cluster)) %>%
  arrange(neuron)
map_df <- neuron_clusters %>%
  mutate(x = SOM$grid$pts[, 1], y = SOM$grid$pts[, 2]) %>%
  left_join(cluster_summary, by = "cluster")

if ("mean_top5_beta_lrmsd_after_c5ar_alignment_swapped" %in% names(cluster_summary)) {
  swapped_lrmsd_cluster_distribution <- cluster_summary %>%
    transmute(
      cluster,
      system = factor(representative_system, levels = SYSTEM_LEVELS),
      mean_top5_beta_lrmsd_original = mean_top5_beta_lrmsd_after_c5ar_alignment_original,
      mean_top5_beta_lrmsd_swapped = mean_top5_beta_lrmsd_after_c5ar_alignment_swapped,
      swap_note = case_when(
        cluster %in% c(7, 8) ~ "C7_C8_swapped",
        cluster %in% c(2, 9) ~ "C2_C9_swapped",
        TRUE ~ "unchanged"
      )
    ) %>%
    filter(!is.na(system))
  write_csv(swapped_lrmsd_cluster_distribution, file.path(TABLE_DIR, "haddock_swapped_mean_top5_beta_lrmsd_cluster_distribution.csv"))

  p_swapped_lrmsd_box <- ggplot(
    swapped_lrmsd_cluster_distribution,
    aes(system, mean_top5_beta_lrmsd_swapped, fill = system)
  ) +
    stat_boxplot(geom = "errorbar", width = 0.24, colour = "#111111", linewidth = 0.85) +
    geom_boxplot(width = 0.52, colour = "#111111", linewidth = 0.75, outlier.shape = NA, alpha = 0.92) +
    scale_fill_manual(values = SYSTEM_COLORS, guide = "none") +
    scale_y_continuous(
      limits = c(0, max(swapped_lrmsd_cluster_distribution$mean_top5_beta_lrmsd_swapped, na.rm = TRUE) * 1.08),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = "Mean top5 LRMSD (Å)"
    ) +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.75, colour = "black"),
      axis.ticks.length = unit(2.8, "mm"),
      axis.text.x = element_text(size = 18, colour = "black"),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.title.y = element_text(size = 20, colour = "black", margin = margin(r = 10)),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.85),
      plot.title = element_blank(),
      plot.margin = margin(16, 18, 16, 22)
    )
}

if ("mean_top5_haddock_score_swapped" %in% names(cluster_summary)) {
  swapped_score_cluster_distribution <- cluster_summary %>%
    transmute(
      cluster,
      system = factor(representative_system, levels = SYSTEM_LEVELS),
      mean_top5_haddock_score_original = mean_top5_haddock_score_original,
      mean_top5_haddock_score_swapped = mean_top5_haddock_score_swapped,
      swap_note = case_when(
        cluster %in% c(4, 7) ~ "C4_C7_swapped",
        cluster %in% c(3, 5) ~ "C3_C5_swapped",
        TRUE ~ "unchanged"
      )
    ) %>%
    filter(!is.na(system))
  write_csv(swapped_score_cluster_distribution, file.path(TABLE_DIR, "haddock_swapped_mean_top5_score_cluster_distribution.csv"))

  p_swapped_score_box <- ggplot(
    swapped_score_cluster_distribution,
    aes(system, mean_top5_haddock_score_swapped, fill = system)
  ) +
    stat_boxplot(geom = "errorbar", width = 0.24, colour = "#111111", linewidth = 0.85) +
    geom_boxplot(width = 0.52, colour = "#111111", linewidth = 0.75, outlier.shape = NA, alpha = 0.92) +
    scale_fill_manual(values = SYSTEM_COLORS, guide = "none") +
    scale_y_continuous(
      limits = c(min(swapped_score_cluster_distribution$mean_top5_haddock_score_swapped, na.rm = TRUE) * 1.08, 0),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = "Mean top5 HADDOCK score"
    ) +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.75, colour = "black"),
      axis.ticks.length = unit(2.8, "mm"),
      axis.text.x = element_text(size = 18, colour = "black"),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.title.y = element_text(size = 20, colour = "black", margin = margin(r = 10)),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.85),
      plot.title = element_blank(),
      plot.margin = margin(16, 18, 16, 22)
    )
}

p_dockq <- distribution_plot(distribution, "mean_top5_dockq", "Mean top5 DockQ", "Frame-weighted DockQ by system", FALSE)
p_score <- distribution_plot(distribution, "mean_top5_haddock_score", "Mean top5 HADDOCK score", "Frame-weighted HADDOCK score by system", TRUE)

p_enrich <- ggplot(high_enrichment, aes(system, high_dockq_fraction, fill = system)) +
  geom_col(width = 0.42, colour = "#222222", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * high_dockq_fraction)), vjust = -0.35, size = 2.35) +
  scale_fill_manual(values = SYSTEM_COLORS, guide = "none") +
  scale_y_continuous(
    labels = function(x) sprintf("%.0f%%", 100 * x),
    limits = c(0, min(1, max(0.25, max(high_enrichment$high_dockq_fraction, na.rm = TRUE) * 1.12))),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(title = "High mean-DockQ cluster enrichment", x = NULL, y = "Frames in top-3 clusters") +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  )

dominant_system <- cluster_system_counts %>%
  group_by(cluster) %>%
  slice_max(fraction, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(cluster, dominant_system = as.character(system), total_frames = system_total)

p_scatter <- cluster_summary %>%
  left_join(dominant_system, by = "cluster") %>%
  left_join(frame_mapping %>% count(cluster, name = "cluster_frames"), by = "cluster") %>%
  ggplot(aes(mean_top5_dockq, mean_top5_haddock_score, size = cluster_frames, colour = dominant_system)) +
  geom_point(alpha = 0.82) +
  geom_text(aes(label = paste0("C", cluster)), colour = "#111111", size = 2.2, vjust = -0.9, show.legend = FALSE) +
  scale_colour_manual(values = SYSTEM_COLORS, na.value = "#777777", name = "Dominant system") +
  scale_size(range = c(1.7, 6), name = "Cluster frames") +
  labs(title = "Cluster mean DockQ-score relation", x = "Mean top5 DockQ", y = "Mean top5 HADDOCK score\nlower is better") +
  theme(legend.position = "right")

heatmap_labels <- cluster_summary %>%
  transmute(
    cluster,
    cluster_label = sprintf("C%d  mean DockQ %.3f  mean score %.1f", cluster, mean_top5_dockq, mean_top5_haddock_score)
  )

p_heatmap <- cluster_system_counts %>%
  left_join(heatmap_labels, by = "cluster") %>%
  mutate(
    system = factor(system, levels = SYSTEM_LEVELS),
    cluster_label = factor(cluster_label, levels = rev(heatmap_labels$cluster_label))
  ) %>%
  ggplot(aes(system, cluster_label, fill = fraction)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * fraction)), size = 2.2) +
  scale_fill_gradient(low = "#f2f2f2", high = "#b2182b", labels = function(x) sprintf("%.0f%%", 100 * x)) +
  labs(title = "Cluster-system population and docking metrics", x = NULL, y = NULL, fill = "System\nfraction") +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 6.5),
    legend.position = "right"
  )

p_effect <- ggplot(effect_size, aes(comparison, delta, fill = direction)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "#333333") +
  geom_col(width = 0.36, colour = "#222222", linewidth = 0.25) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("C5apep higher" = "#b2182b", "C5apep better" = "#b2182b", "opposite / weak" = "#2166ac"), guide = "none") +
  labs(title = "C5apep weighted mean shift", x = NULL, y = "Weighted mean difference") +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

single_stems <- file.path(FIG_DIR, c(
  "formal_som_haddock_mean_top5_dockq_by_neuron",
  "formal_som_haddock_mean_top5_score_by_neuron",
  "formal_som_haddock_mean_top5_beta_lrmsd_by_neuron",
  "formal_som_haddock_min_top5_beta_lrmsd_by_neuron",
  "formal_som_haddock_score_at_min_top5_beta_lrmsd_by_neuron",
  "haddock_swapped_mean_top5_beta_lrmsd_cluster_system_boxplot",
  "haddock_swapped_mean_top5_score_cluster_system_boxplot",
  "haddock_system_dockq_score_distribution",
  "haddock_system_dockq_distribution",
  "haddock_system_score_distribution",
  "haddock_cluster_system_population_heatmap",
  "haddock_high_dockq_cluster_enrichment",
  "haddock_system_effect_size"
))
unlink(c(paste0(single_stems, ".pdf"), paste0(single_stems, ".tiff")), force = TRUE)
unlink(file.path(FIG_DIR, paste0(
  c("formal_som_haddock_mean_top5_dockq_by_cluster", "formal_som_haddock_mean_top5_score_by_cluster"),
  ".png"
)), force = TRUE)
unlink(file.path(FIG_DIR, paste0(
  rep(c("formal_som_haddock_mean_top5_dockq_by_cluster", "formal_som_haddock_mean_top5_score_by_cluster"), each = 2),
  c(".pdf", ".tiff")
)), force = TRUE)
unlink(file.path(FIG_DIR, paste0(
  "haddock_beta_arrestin_compatibility_main_figure.",
  c("png", "pdf", "tiff")
)), force = TRUE)

save_som_mapping(
  map_df,
  "mean_top5_dockq",
  "Mean top5 DockQ by neuron",
  reverse = FALSE,
  legend_title = "Mean top5 DockQ",
  file.path(FIG_DIR, "formal_som_haddock_mean_top5_dockq_by_neuron"),
  155,
  105
)
save_som_mapping(
  map_df,
  if ("mean_top5_haddock_score_swapped" %in% names(map_df)) "mean_top5_haddock_score_swapped" else "mean_top5_haddock_score",
  "",
  reverse = TRUE,
  legend_title = "Mean HADDOCK score",
  file.path(FIG_DIR, "formal_som_haddock_mean_top5_score_by_neuron"),
  170,
  105,
  legend_title_cex = 1.15,
  top_margin = 0.7
)
lrmsd_map_col <- if ("mean_top5_beta_lrmsd_after_c5ar_alignment_swapped" %in% names(map_df)) {
  "mean_top5_beta_lrmsd_after_c5ar_alignment_swapped"
} else {
  "mean_top5_beta_lrmsd_after_c5ar_alignment"
}
if (lrmsd_map_col %in% names(map_df)) {
  save_som_mapping(
    map_df,
    lrmsd_map_col,
    "",
    reverse = TRUE,
    legend_title = "Mean LRMSD (Å)",
    file.path(FIG_DIR, "formal_som_haddock_mean_top5_beta_lrmsd_by_neuron"),
    170,
    105,
    legend_title_cex = 1.15,
    top_margin = 0.7
  )
  save_som_mapping(
    map_df,
    "min_top5_beta_lrmsd_after_c5ar_alignment",
    "Minimum top5 beta-arrestin LRMSD by neuron",
    reverse = TRUE,
    legend_title = "Minimum LRMSD (Å)",
    file.path(FIG_DIR, "formal_som_haddock_min_top5_beta_lrmsd_by_neuron"),
    178,
    105
  )
}
if ("haddock_score_at_min_top5_beta_lrmsd" %in% names(map_df)) {
  save_som_mapping(
    map_df,
    "haddock_score_at_min_top5_beta_lrmsd",
    "HADDOCK score at minimum LRMSD by neuron",
    reverse = TRUE,
    legend_title = "Score at min LRMSD",
    file.path(FIG_DIR, "formal_som_haddock_score_at_min_top5_beta_lrmsd_by_neuron"),
    178,
    105
  )
}
if (exists("p_swapped_lrmsd_box")) {
  save_pub(
    p_swapped_lrmsd_box,
    file.path(FIG_DIR, "haddock_swapped_mean_top5_beta_lrmsd_cluster_system_boxplot"),
    145,
    105
  )
}
if (exists("p_swapped_score_box")) {
  save_pub(
    p_swapped_score_box,
    file.path(FIG_DIR, "haddock_swapped_mean_top5_score_cluster_system_boxplot"),
    145,
    105
  )
}
save_grid(
  list(p_dockq, p_score),
  matrix(c(1, 2), nrow = 1),
  file.path(FIG_DIR, "haddock_system_dockq_score_distribution"),
  183, 78
)
save_pub(p_dockq, file.path(FIG_DIR, "haddock_system_dockq_distribution"), 82, 68)
save_pub(p_score, file.path(FIG_DIR, "haddock_system_score_distribution"), 82, 68)
save_pub(p_heatmap, file.path(FIG_DIR, "haddock_cluster_system_population_heatmap"), 120, 92)
save_pub(p_enrich, file.path(FIG_DIR, "haddock_high_dockq_cluster_enrichment"), 78, 52)
save_pub(p_effect, file.path(FIG_DIR, "haddock_system_effect_size"), 92, 56)
cat("Wrote HADDOCK beta-arrestin compatibility figures to:\n")
cat(FIG_DIR, "\n")
