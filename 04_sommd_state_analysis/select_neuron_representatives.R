#!/usr/bin/env Rscript

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    root = "/mnt/e/work/modeling/c5ar/md",
    profile = "formal",
    out_csv = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--root", "--profile", "--out-csv")) {
      if (i == length(args)) stop("Missing value for ", key)
      value <- args[[i + 1]]
      if (key == "--root") out$root <- value
      if (key == "--profile") out$profile <- value
      if (key == "--out-csv") out$out_csv <- value
      i <- i + 2
    } else {
      stop("Unknown argument: ", key)
    }
  }
  out
}

args <- parse_args()
root <- normalizePath(args$root, mustWork = TRUE)
sommd_root <- file.path(root, "article/analysis/sommd")
run_dir <- file.path(sommd_root, "runs", args$profile, "activation_core")
table_dir <- file.path(sommd_root, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- args$out_csv
if (is.null(out_csv)) {
  out_csv <- file.path(table_dir, paste0(args$profile, "_neuron_representative_frames.csv"))
}

som_path <- file.path(run_dir, "SOM.rds")
training_path <- file.path(run_dir, "feature_matrix_training.rds")
frame_map_path <- file.path(run_dir, "frame_mapping.csv")
neuron_cluster_path <- file.path(run_dir, "neuron_clusters.csv")

for (path in c(som_path, training_path, frame_map_path, neuron_cluster_path)) {
  if (!file.exists(path)) stop("Missing required file: ", path)
}

SOM <- readRDS(som_path)
X <- readRDS(training_path)
frame_map <- read.csv(frame_map_path, stringsAsFactors = FALSE)
neuron_clusters <- read.csv(neuron_cluster_path, stringsAsFactors = FALSE)

if (nrow(X) != nrow(frame_map)) {
  stop("Training feature matrix and frame_mapping row counts differ")
}
if (!"neuron" %in% colnames(frame_map)) {
  stop("frame_mapping.csv does not contain a neuron column")
}

codes <- SOM$codes[[1]]
n_neurons <- nrow(codes)
systems <- c("apo", "BM213", "C5apep")
cluster_by_neuron <- setNames(neuron_clusters$cluster, neuron_clusters$neuron)

rows <- vector("list", n_neurons)
for (neuron_id in seq_len(n_neurons)) {
  frames <- which(frame_map$neuron == neuron_id)
  system_counts <- table(factor(frame_map$system[frames], levels = systems))
  frame_count <- length(frames)
  dominant_system <- NA_character_
  dominant_fraction <- NA_real_
  representative_som_frame <- NA_integer_
  representative_frame_neuron <- NA_integer_
  representative_frame_cluster <- NA_integer_
  representative_assigned_to_neuron <- NA
  selection_scope <- "empty"
  distance_to_neuron <- NA_real_
  meta <- as.list(setNames(rep(NA_character_, 6), c(
    "system", "replicate", "pdb", "source_trajectory",
    "downsampled_trajectory", "profile"
  )))
  local_frame_index <- NA_integer_
  source_frame_index_estimate <- NA_integer_
  time_ps <- NA_real_
  time_ns <- NA_real_
  global_frame <- NA_integer_

  if (frame_count > 0) {
    dominant_system <- names(system_counts)[which.max(system_counts)]
    dominant_fraction <- as.numeric(max(system_counts)) / frame_count
    candidate_frames <- frames
    selection_scope <- "assigned_neuron_nearest_codebook"
    representative_assigned_to_neuron <- TRUE
  } else {
    candidate_frames <- seq_len(nrow(X))
    selection_scope <- "global_nearest_codebook_empty_neuron"
    representative_assigned_to_neuron <- FALSE
  }

  if (length(candidate_frames) > 0) {
    code_matrix <- matrix(codes[neuron_id, ], nrow = length(candidate_frames), ncol = ncol(codes), byrow = TRUE)
    distances <- sqrt(rowSums((X[candidate_frames, , drop = FALSE] - code_matrix) ^ 2))
    best_row <- candidate_frames[[which.min(distances)]]
    representative_som_frame <- frame_map$som_frame[[best_row]]
    representative_frame_neuron <- frame_map$neuron[[best_row]]
    representative_frame_cluster <- frame_map$cluster[[best_row]]
    distance_to_neuron <- min(distances)
    meta <- as.list(frame_map[best_row, c(
      "system", "replicate", "pdb", "source_trajectory",
      "downsampled_trajectory", "profile"
    )])
    local_frame_index <- frame_map$local_frame_index[[best_row]]
    source_frame_index_estimate <- frame_map$source_frame_index_estimate[[best_row]]
    time_ps <- frame_map$time_ps[[best_row]]
    time_ns <- frame_map$time_ns[[best_row]]
    global_frame <- frame_map$global_frame[[best_row]]
  }

  rows[[neuron_id]] <- data.frame(
    neuron = neuron_id,
    cluster = as.integer(cluster_by_neuron[[as.character(neuron_id)]]),
    frame_count = frame_count,
    apo_frames = as.integer(system_counts[["apo"]]),
    BM213_frames = as.integer(system_counts[["BM213"]]),
    C5apep_frames = as.integer(system_counts[["C5apep"]]),
    dominant_system = dominant_system,
    dominant_system_fraction = dominant_fraction,
    selection_scope = selection_scope,
    representative_assigned_to_neuron = representative_assigned_to_neuron,
    representative_som_frame = representative_som_frame,
    representative_frame_neuron = representative_frame_neuron,
    representative_frame_cluster = representative_frame_cluster,
    distance_to_neuron = distance_to_neuron,
    global_frame = global_frame,
    profile = meta$profile,
    system = meta$system,
    replicate = meta$replicate,
    local_frame_index = local_frame_index,
    source_frame_index_estimate = source_frame_index_estimate,
    time_ps = time_ps,
    time_ns = time_ns,
    pdb = meta$pdb,
    source_trajectory = meta$source_trajectory,
    downsampled_trajectory = meta$downsampled_trajectory,
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE)
writeLines(out_csv)
