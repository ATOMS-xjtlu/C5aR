#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SOMMD)
  library(kohonen)
})

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    token <- sub("^--", "", token)
    if (grepl("=", token, fixed = TRUE)) {
      parts <- strsplit(token, "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      value <- paste(parts[-1], collapse = "=")
    } else {
      key <- token
      if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
        value <- "true"
      } else {
        i <- i + 1
        value <- argv[[i]]
      }
    }
    key <- gsub("-", "_", key)
    args[[key]] <- value
    i <- i + 1
  }
  args
}

arg <- parse_args(commandArgs(trailingOnly = TRUE))

if (!is.null(arg[["help"]])) {
  cat(paste(
    "Usage: run_activation_core_som.R",
    "--feature-matrix feature_matrix.csv",
    "--frame-metadata frame_metadata.csv",
    "--feature-metadata feature_metadata.csv",
    "--out outdir",
    "[--grid 8x8] [--clusters 8] [--rlen 500] [--seed 20260603]",
    "[--toroidal true] [--scale true]",
    "\n"
  ))
  quit(status = 0)
}

get_arg <- function(name, default = NULL, required = FALSE) {
  if (!is.null(arg[[name]]) && nzchar(arg[[name]])) return(arg[[name]])
  if (required) stop("Missing required argument --", gsub("_", "-", name))
  default
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

split_grid <- function(x) {
  x <- tolower(gsub(" ", "", x))
  parts <- strsplit(x, "x", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop("--grid must look like 8x8")
  dims <- as.integer(parts)
  if (any(is.na(dims)) || any(dims < 2)) stop("Invalid --grid: ", x)
  dims
}

write_plot <- function(filename, expr, width = 1600, height = 1200) {
  png(filename, width = width, height = height, res = 200)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

feature_matrix_path <- normalizePath(get_arg("feature_matrix", required = TRUE), mustWork = TRUE)
frame_metadata_path <- normalizePath(get_arg("frame_metadata", required = TRUE), mustWork = TRUE)
feature_metadata_path <- normalizePath(get_arg("feature_metadata", required = TRUE), mustWork = TRUE)
outdir <- get_arg("out", "sommd_activation_core_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cache_dir <- file.path(outdir, ".cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(XDG_CACHE_HOME = normalizePath(cache_dir, mustWork = FALSE))

grid_dims <- split_grid(get_arg("grid", "8x8"))
toroidal <- as_bool(get_arg("toroidal", "true"))
clusters <- as.integer(get_arg("clusters", "8"))
rlen <- as.integer(get_arg("rlen", "500"))
seed <- as.integer(get_arg("seed", "20260603"))
scale_features <- as_bool(get_arg("scale", "true"))

if (is.na(clusters) || clusters < 2) stop("--clusters must be >= 2")
if (is.na(rlen) || rlen < 1) stop("--rlen must be positive")

feature_df <- read.csv(feature_matrix_path, check.names = FALSE)
metadata <- read.csv(frame_metadata_path, check.names = FALSE)
feature_metadata <- read.csv(feature_metadata_path, check.names = FALSE)
X_raw <- data.matrix(feature_df)

if (nrow(X_raw) < 2) stop("Feature matrix has fewer than 2 frames")
if (nrow(metadata) != nrow(X_raw)) stop("frame metadata rows do not match feature matrix rows")
if (ncol(X_raw) != 91) stop("Expected 91 activation-core distance features, found ", ncol(X_raw))
if (any(!is.finite(X_raw))) stop("Feature matrix contains non-finite values")

center <- rep(0, ncol(X_raw))
scale <- rep(1, ncol(X_raw))
X <- X_raw
if (scale_features) {
  center <- colMeans(X_raw)
  scale <- apply(X_raw, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  X <- sweep(sweep(X_raw, 2, center, "-"), 2, scale, "/")
}

message("Feature matrix: ", nrow(X), " frames x ", ncol(X), " features")
set.seed(seed)
SOM <- kohonen::som(
  X,
  grid = kohonen::somgrid(
    grid_dims[[1]], grid_dims[[2]], "hexagonal",
    neighbourhood.fct = "gaussian",
    toroidal = toroidal
  ),
  dist.fcts = "euclidean",
  rlen = rlen,
  mode = "pbatch"
)

saveRDS(SOM, file.path(outdir, "SOM.rds"))
saveRDS(
  list(
    feature = "activation_core_ca_distances",
    feature_matrix = feature_matrix_path,
    frame_metadata = frame_metadata_path,
    feature_metadata = feature_metadata_path,
    grid = grid_dims,
    toroidal = toroidal,
    clusters = clusters,
    rlen = rlen,
    seed = seed,
    scaled = scale_features,
    center = center,
    scale = scale
  ),
  file.path(outdir, "feature_info.rds")
)
saveRDS(X_raw, file.path(outdir, "feature_matrix_raw.rds"))
saveRDS(X, file.path(outdir, "feature_matrix_training.rds"))
write.csv(
  data.frame(feature_name = colnames(X_raw), center = center, scale = scale),
  file.path(outdir, "feature_scaling.csv"),
  row.names = FALSE
)

write_plot(file.path(outdir, "01_umatrix.png"), {
  plot(SOM, type = "dist.neighbours", heatkey = TRUE, shape = "straight", main = "U-Matrix")
})

max_clusters <- min(30, max(3, nrow(SOM$codes[[1]]) - 1))
sil_score <- silhouette.score(SOM, clust_method = "complete", interval = seq(2, max_clusters))
write.csv(sil_score, file.path(outdir, "02_silhouette_score.csv"), row.names = FALSE)
write_plot(file.path(outdir, "02_silhouette_score.png"), {
  plot(sil_score, type = "b", pch = 19, lwd = 1, xlab = "Number of clusters", ylab = "Average silhouettes")
  abline(v = clusters, col = "red", lty = 3, lwd = 2)
})

SOM_hc <- cutree(hclust(dist(SOM$codes[[1]], method = "euclidean"), method = "complete"), clusters)
colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ffff88", "#6a3d9a",
            "#a0451f", "#96c3dc", "#fbb25c", "#ff7f00", "#bea0cc",
            "#747474", "#f88587", "#a4db77", "#b15928", "#8dd3c7")
if (clusters > length(colors)) {
  colors <- grDevices::colorRampPalette(colors)(clusters)
}

write.csv(
  data.frame(neuron = seq_along(SOM_hc), cluster = SOM_hc),
  file.path(outdir, "neuron_clusters.csv"),
  row.names = FALSE
)

write_plot(file.path(outdir, "03_som_clusters.png"), {
  plot(SOM, type = "mapping", bgcol = colors[SOM_hc], col = rgb(0, 0, 0, 0), shape = "straight", main = "")
  kohonen::add.cluster.boundaries(SOM, SOM_hc, lwd = 4)
  som.add.numbers(SOM, scale = 0.45, col = "black")
})

som_frame <- seq_len(nrow(X))
neuron <- SOM$unit.classif
cluster <- SOM_hc[neuron]
frame_map <- cbind(metadata, data.frame(som_frame = som_frame, neuron = neuron, cluster = cluster))
write.csv(frame_map, file.path(outdir, "frame_mapping.csv"), row.names = FALSE)
write.csv(
  as.data.frame(table(cluster = frame_map$cluster)),
  file.path(outdir, "cluster_populations.csv"),
  row.names = FALSE
)

representatives <- data.frame()
codes <- SOM$codes[[1]]
for (cl in seq_len(clusters)) {
  frames <- which(cluster == cl)
  if (length(frames) == 0) next
  unit <- neuron[frames]
  d <- vapply(seq_along(frames), function(i) {
    sqrt(sum((X[frames[[i]], ] - codes[unit[[i]], ]) ^ 2))
  }, numeric(1))
  best <- frames[[which.min(d)]]
  representatives <- rbind(
    representatives,
    data.frame(
      cluster = cl,
      representative_neuron = neuron[[best]],
      representative_frame = best,
      distance_to_neuron = min(d)
    )
  )
}
write.csv(representatives, file.path(outdir, "representative_frames.csv"), row.names = FALSE)

transition <- matrix(0L, nrow = clusters, ncol = clusters)
rownames(transition) <- paste0("from_", seq_len(clusters))
colnames(transition) <- paste0("to_", seq_len(clusters))
groups <- interaction(frame_map$system, frame_map$replicate, drop = TRUE)
for (g in levels(groups)) {
  seg <- which(groups == g)
  if (length(seg) < 2) next
  a <- cluster[seg[-length(seg)]]
  b <- cluster[seg[-1]]
  keep <- !is.na(a) & !is.na(b)
  for (j in which(keep)) transition[a[[j]], b[[j]]] <- transition[a[[j]], b[[j]]] + 1L
}
write.csv(transition, file.path(outdir, "transition_matrix.csv"))

run_info <- c(
  "Arguments:",
  paste(capture.output(str(arg)), collapse = "\n"),
  "",
  "Resolved inputs:",
  paste("feature_matrix:", feature_matrix_path),
  paste("frame_metadata:", frame_metadata_path),
  paste("feature_metadata:", feature_metadata_path),
  paste("training_matrix:", nrow(X), "frames x", ncol(X), "features"),
  paste("scaled:", scale_features),
  "",
  "sessionInfo:",
  capture.output(sessionInfo())
)
writeLines(run_info, file.path(outdir, "run_parameters.txt"))

message("Activation-core SOM complete: ", normalizePath(outdir, mustWork = FALSE))
