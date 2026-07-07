#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(kohonen)
})

ROOT <- "/mnt/e/work/modeling/c5ar/md"
SOMMD_ROOT <- file.path(ROOT, "article/analysis/sommd")
RUN_DIR <- file.path(SOMMD_ROOT, "runs/formal/activation_core")
FIG_DIR <- file.path(SOMMD_ROOT, "figures")
TABLE_DIR <- file.path(SOMMD_ROOT, "tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SOM <- readRDS(file.path(RUN_DIR, "SOM.rds"))
neuron_clusters <- read.csv(file.path(RUN_DIR, "neuron_clusters.csv"))
helix <- read.csv(file.path(TABLE_DIR, "formal_loop_helix_by_neuron.csv"))

segments <- c("ICL2", "ICL3", "H8")
palette <- grDevices::colorRampPalette(c("#f7fbff", "#9ecae1", "#3182bd", "#08306b"))(101)
na_col <- "#e5e5e5"

value_to_color <- function(values) {
  out <- rep(na_col, length(values))
  ok <- is.finite(values)
  idx <- pmax(1, pmin(101, floor(values[ok] * 100) + 1))
  out[ok] <- palette[idx]
  out
}

draw_colorbar <- function(
  label = "Mean helix fraction",
  label_cex = 0.95,
  axis_cex = 0.9,
  draw_box = TRUE,
  label_side = c("top", "right")
) {
  label_side <- match.arg(label_side)
  par(mar = c(0.8, 0.1, 2.0, 0.2))
  plot.new()
  plot.window(xlim = c(0, 2.75), ylim = c(-0.05, 1.05), xaxs = "i", yaxs = "i")
  y <- seq(0, 1, length.out = 101)
  bar_left <- 0.22
  bar_right <- 0.62
  axis_x <- bar_right
  tick_x <- 0.76
  label_x <- 0.90
  title_x <- 1.55
  for (i in seq_len(100)) {
    rect(bar_left, y[i], bar_right, y[i + 1], col = palette[i], border = NA)
  }
  rect(bar_left, 0, bar_right, 1, border = "#222222", lwd = 0.9)
  segments(axis_x, 0, axis_x, 1, col = "#111111", lwd = 1.0)
  ticks <- seq(0, 1, 0.25)
  segments(axis_x, ticks, tick_x, ticks, col = "#111111", lwd = 1.0)
  text(label_x, ticks, labels = sprintf("%.2f", ticks), adj = c(0, 0.5), cex = axis_cex)
  if (label_side == "right") {
    text(title_x, 0.5, labels = label, srt = 270, adj = c(0.5, 0.5), cex = label_cex, xpd = NA)
  } else {
    text(0.16, 1.035, labels = label, adj = c(0, 0), cex = label_cex, xpd = NA)
  }
  if (draw_box) {
    rect(0.12, -0.035, 1.50, 1.025, border = "#111111", lwd = 1.0)
  }
}

plot_segment <- function(segment, main = NULL) {
  sub <- helix[helix$segment == segment, ]
  sub <- sub[match(neuron_clusters$neuron, sub$neuron), ]
  colors <- value_to_color(sub$mean_helix_fraction)
  plot(
    SOM,
    type = "mapping",
    bgcol = colors,
    col = rgb(0, 0, 0, 0),
    shape = "straight",
    main = ifelse(is.null(main), paste0(segment, " helix content"), main)
  )
  kohonen::add.cluster.boundaries(SOM, neuron_clusters$cluster, lwd = 4)
}

for (segment in segments) {
  outfile <- file.path(FIG_DIR, paste0("formal_som_", segment, "_helix_by_neuron.png"))
  png(outfile, width = 2200, height = 1450, res = 220)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1.75))
  par(mar = c(0.5, 0.5, 3.1, 0.2), cex.main = 1.65)
  plot_segment(segment)
  draw_colorbar(label_cex = 1.45, axis_cex = 1.18, draw_box = FALSE, label_side = "right")
  dev.off()
}

png(file.path(FIG_DIR, "formal_som_loop_helix_by_neuron.png"), width = 3650, height = 1250, res = 220)
layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(1, 1, 1, 0.62))
for (segment in segments) {
  par(mar = c(0.4, 0.4, 2.8, 0.2), cex.main = 1.28)
  plot_segment(segment)
}
draw_colorbar(label_cex = 1.22, axis_cex = 0.9, draw_box = FALSE)
dev.off()

writeLines(c(
  file.path(FIG_DIR, "formal_som_ICL2_helix_by_neuron.png"),
  file.path(FIG_DIR, "formal_som_ICL3_helix_by_neuron.png"),
  file.path(FIG_DIR, "formal_som_H8_helix_by_neuron.png"),
  file.path(FIG_DIR, "formal_som_loop_helix_by_neuron.png")
))
