#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SOMMD)
  library(kohonen)
})

ROOT <- "/mnt/e/work/modeling/c5ar/md"
SOMMD_ROOT <- file.path(ROOT, "article/analysis/sommd")
RUN_DIR <- file.path(SOMMD_ROOT, "runs/formal/activation_core")
FIG_DIR <- file.path(SOMMD_ROOT, "figures")
TABLE_DIR <- file.path(SOMMD_ROOT, "tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

SOM <- readRDS(file.path(RUN_DIR, "SOM.rds"))
neuron_clusters <- read.csv(file.path(RUN_DIR, "neuron_clusters.csv"))
frame_map <- read.csv(file.path(RUN_DIR, "frame_mapping.csv"))

cluster_ids <- sort(unique(neuron_clusters$cluster))
system_order <- c("apo", "BM213", "C5apep")
cluster_cols <- c(
  "#1f78b4", "#33a02c", "#e31a1c", "#ffff88", "#6a3d9a",
  "#a0451f", "#96c3dc", "#fbb25c", "#ff7f00", "#bea0cc",
  "#747474", "#f88587", "#a4db77", "#b15928", "#8dd3c7"
)
cluster_cols <- cluster_cols[seq_along(cluster_ids)]
names(cluster_cols) <- as.character(cluster_ids)
system_cols <- c(apo = "#4d4d4d", BM213 = "#00BFC4", C5apep = "#e68613")

add_cluster_labels <- function(cex = 0.58, offset = c(0, 0), font = 2, include_neuron = FALSE) {
  xy <- SOM$grid$pts
  labels <- paste0("C", neuron_clusters$cluster)
  if (include_neuron) {
    labels <- paste0(neuron_clusters$neuron, "\n", labels)
  }
  text(
    xy[, 1] + offset[1],
    xy[, 2] + offset[2],
    labels = labels,
    cex = cex,
    font = font,
    col = "#111111"
  )
}

plot_cluster_background <- function(
  main = "",
  label_clusters = TRUE,
  label_cex = 0.58,
  label_offset = c(0, 0),
  include_neuron = FALSE
) {
  bg <- cluster_cols[as.character(neuron_clusters$cluster)]
  plot(
    SOM,
    type = "mapping",
    bgcol = bg,
    col = rgb(0, 0, 0, 0),
    shape = "straight",
    main = main
  )
  kohonen::add.cluster.boundaries(SOM, neuron_clusters$cluster, lwd = 4)
  if (label_clusters) {
    add_cluster_labels(cex = label_cex, offset = label_offset, include_neuron = include_neuron)
  }
}

legend_labels <- paste0("Cluster ", cluster_ids)

png(file.path(FIG_DIR, "formal_som_clusters_with_legend.png"), width = 1900, height = 1780, res = 220)
layout(matrix(c(1, 2), nrow = 2), heights = c(5.1, 1.35))
par(mar = c(0.2, 0.5, 2.0, 0.5))
plot_cluster_background(
  "formal SOM clusters",
  label_clusters = TRUE,
  label_cex = 0.64,
  include_neuron = TRUE
)
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = legend_labels,
  fill = cluster_cols,
  border = "#333333",
  bty = "n",
  ncol = 3,
  cex = 1.52
)
dev.off()
file.copy(
  file.path(FIG_DIR, "formal_som_clusters_with_legend.png"),
  file.path(RUN_DIR, "03_som_clusters_with_legend.png"),
  overwrite = TRUE
)

pop <- as.data.frame(table(frame_map$system, frame_map$cluster), stringsAsFactors = FALSE)
colnames(pop) <- c("system", "cluster", "frames")
pop$cluster <- as.integer(as.character(pop$cluster))
pop$system <- factor(pop$system, levels = system_order)
pop <- pop[order(pop$system, pop$cluster), ]
totals <- tapply(pop$frames, pop$system, sum)
pop$percentage <- 100 * pop$frames / totals[as.character(pop$system)]
write.csv(pop, file.path(TABLE_DIR, "formal_cluster_population_percent_by_system.csv"), row.names = FALSE)

bar_mat <- matrix(
  0,
  nrow = length(system_order),
  ncol = length(cluster_ids),
  dimnames = list(system_order, paste("Cluster", cluster_ids))
)
for (i in seq_len(nrow(pop))) {
  bar_mat[as.character(pop$system[i]), paste("Cluster", pop$cluster[i])] <- pop$percentage[i]
}
plot_bar_mat <- bar_mat
plot_bar_mat[plot_bar_mat == 0] <- NA

draw_cluster_population_bars <- function(
  values,
  main,
  cex_axis = 1.18,
  cex_lab = 1.28,
  cex_main = 1.18,
  cex_names = 1.16,
  legend_cex = 1.55,
  legend_inset = c(-0.12, 0.02)
) {
  n_systems <- nrow(values)
  n_clusters <- ncol(values)
  cluster_centers <- seq_len(n_clusters) * 1.35
  bar_width <- 0.26
  offsets <- (seq_len(n_systems) - (n_systems + 1) / 2) * 0.33
  x_pos <- matrix(rep(cluster_centers, each = n_systems), nrow = n_systems) +
    matrix(rep(offsets, times = n_clusters), nrow = n_systems)
  x_limits <- range(cluster_centers) + c(-0.85, 0.85)

  plot.new()
  plot.window(xlim = x_limits, ylim = c(0, 100), xaxs = "i", yaxs = "i")
  title(main = main, cex.main = cex_main)

  for (j in seq_len(n_clusters)) {
    for (i in seq_len(n_systems)) {
      height <- values[i, j]
      if (is.finite(height) && height > 0) {
        rect(
          x_pos[i, j] - bar_width / 2,
          0,
          x_pos[i, j] + bar_width / 2,
          height,
          col = system_cols[rownames(values)[i]],
          border = "black",
          lwd = 0.9
        )
      }
    }
  }

  for (j in seq_len(n_clusters - 1)) {
    x_sep <- (cluster_centers[j] + cluster_centers[j + 1]) / 2
    segments(x_sep, 0, x_sep, 100, lty = 2, col = "#666666", lwd = 0.9)
  }

  axis(2, las = 1, lwd = 1.2, cex.axis = cex_axis)
  segments(x_limits[1], 0, x_limits[2], 0, lwd = 1.2, col = "black", xpd = FALSE)
  mtext("Population (%)", side = 2, line = 3.7, cex = cex_lab)
  text(
    cluster_centers,
    -6.5,
    labels = paste0("C", cluster_ids),
    cex = cex_names,
    xpd = NA
  )

  par(xpd = NA)
  legend(
    "topright",
    inset = legend_inset,
    legend = system_order,
    fill = system_cols[system_order],
    border = "black",
    bty = "n",
    cex = legend_cex
  )
  par(xpd = FALSE)
  invisible(x_pos)
}

png(file.path(FIG_DIR, "formal_cluster_population_barplot.png"), width = 3000, height = 1400, res = 220)
par(mar = c(6.0, 5.8, 2.6, 10.8), xpd = FALSE)
draw_cluster_population_bars(
  bar_mat,
  main = "Formal cluster population by system",
  cex_axis = 1.22,
  cex_lab = 1.34,
  cex_main = 1.22,
  cex_names = 1.20,
  legend_cex = 1.55,
  legend_inset = c(-0.12, 0.02)
)
dev.off()

neuron_grid <- expand.grid(system = system_order, neuron = seq_len(nrow(SOM$grid$pts)))
neuron_counts <- as.data.frame(table(frame_map$system, frame_map$neuron), stringsAsFactors = FALSE)
colnames(neuron_counts) <- c("system", "neuron", "frames")
neuron_counts$neuron <- as.integer(as.character(neuron_counts$neuron))
neuron_occupancy <- merge(neuron_grid, neuron_counts, by = c("system", "neuron"), all.x = TRUE)
neuron_occupancy$frames[is.na(neuron_occupancy$frames)] <- 0
neuron_occupancy <- merge(neuron_occupancy, neuron_clusters, by = "neuron", all.x = TRUE)
neuron_occupancy <- neuron_occupancy[order(factor(neuron_occupancy$system, levels = system_order), neuron_occupancy$neuron), ]
write.csv(neuron_occupancy, file.path(TABLE_DIR, "formal_neuron_occupancy_by_system.csv"), row.names = FALSE)

add_occupancy_points <- function(system) {
  sub <- neuron_occupancy[neuron_occupancy$system == system & neuron_occupancy$frames > 0, ]
  if (nrow(sub) == 0) {
    return(invisible(NULL))
  }
  xy <- SOM$grid$pts[sub$neuron, , drop = FALSE]
  max_frames <- max(neuron_occupancy$frames)
  point_cex <- 0.28 + 1.62 * sqrt(sub$frames / max_frames)
  points(
    xy[, 1],
    xy[, 2],
    pch = 21,
    bg = "#ff00cc",
    col = "#8b005f",
    lwd = 0.7,
    cex = point_cex
  )
}

png(file.path(FIG_DIR, "formal_som_neuron_occupancy_by_system.png"), width = 2400, height = 1000, res = 220)
par(mfrow = c(1, 3), mar = c(0.2, 0.3, 2.0, 0.3), oma = c(0.2, 0, 2.0, 0))
for (system in system_order) {
  plot_cluster_background(system, label_clusters = FALSE)
  add_occupancy_points(system)
  add_cluster_labels(cex = 0.42, offset = c(-0.24, 0.26), font = 2)
}
mtext("Magenta circle area is proportional to neuron frame count", outer = TRUE, side = 3, line = 0.1, cex = 0.95)
dev.off()

for (system in system_order) {
  outfile <- file.path(FIG_DIR, paste0("formal_som_neuron_occupancy_", system, ".png"))
  png(outfile, width = 1500, height = 1300, res = 220)
  par(mar = c(0.5, 0.5, 2.0, 0.5))
  plot_cluster_background(system, label_clusters = FALSE)
  add_occupancy_points(system)
  add_cluster_labels(cex = 0.50, offset = c(-0.24, 0.26), font = 2)
  dev.off()
}

png(file.path(FIG_DIR, "formal_som_population_composite.png"), width = 2600, height = 1800, res = 220)
layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = TRUE), heights = c(1.15, 1.0))
par(mar = c(4.8, 5.0, 2.0, 4.4), xpd = FALSE)
draw_cluster_population_bars(
  bar_mat,
  main = "Formal cluster population by system",
  cex_axis = 1.0,
  cex_lab = 1.1,
  cex_main = 1.05,
  cex_names = 0.95,
  legend_cex = 0.98,
  legend_inset = c(0.01, 0.03)
)
for (system in system_order) {
  par(mar = c(0.2, 0.2, 2.0, 0.2))
  plot_cluster_background(system, label_clusters = FALSE)
  add_occupancy_points(system)
  add_cluster_labels(cex = 0.38, offset = c(-0.24, 0.26), font = 2)
}
mtext("Magenta circle area is proportional to neuron frame count", outer = TRUE, side = 1, line = -0.1, cex = 0.8)
dev.off()

writeLines(c(
  file.path(FIG_DIR, "formal_som_clusters_with_legend.png"),
  file.path(RUN_DIR, "03_som_clusters_with_legend.png"),
  file.path(FIG_DIR, "formal_cluster_population_barplot.png"),
  file.path(FIG_DIR, "formal_som_neuron_occupancy_by_system.png"),
  file.path(FIG_DIR, "formal_som_population_composite.png"),
  file.path(TABLE_DIR, "formal_cluster_population_percent_by_system.csv"),
  file.path(TABLE_DIR, "formal_neuron_occupancy_by_system.csv")
))
