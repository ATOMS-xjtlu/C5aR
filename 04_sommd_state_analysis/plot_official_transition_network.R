#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SOMMD)
  library(igraph)
})

ROOT <- "/mnt/e/work/modeling/c5ar/md"
SOMMD_ROOT <- file.path(ROOT, "article/analysis/sommd")
SYSTEM_ORDER <- c("apo", "BM213", "C5apep")
COL_SCALE <- c(
  "#1f78b4", "#33a02c", "#e31a1c", "#ffff88", "#6a3d9a",
  "#a0451f", "#96c3dc", "#fbb25c", "#ff7f00", "#bea0cc",
  "#747474", "#f88587", "#a4db77", "#8dd3c7", "#fb8072"
)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    profile = "formal",
    run_dir = NULL,
    out_figures = file.path(SOMMD_ROOT, "figures"),
    out_tables = file.path(SOMMD_ROOT, "tables"),
    seed = 20260603
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(paste("Unexpected argument:", key))
    }
    name <- substring(key, 3)
    if (!(name %in% names(out))) {
      stop(paste("Unknown argument:", key))
    }
    if (i == length(args)) {
      stop(paste("Missing value for", key))
    }
    out[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  out$seed <- as.integer(out$seed)
  if (is.null(out$run_dir)) {
    out$run_dir <- file.path(SOMMD_ROOT, "runs", out$profile, "activation_core")
  }
  out
}

make_transition_matrix <- function(SOM, frame_map, systems = NULL) {
  n_neurons <- nrow(SOM$codes[[1]])
  out <- matrix(0, nrow = n_neurons, ncol = n_neurons)
  rownames(out) <- paste("N_", seq_len(n_neurons), sep = "")
  colnames(out) <- paste("N_", seq_len(n_neurons), sep = "")
  sub <- frame_map
  if (!is.null(systems)) {
    sub <- sub[sub$system %in% systems, , drop = FALSE]
  }
  groups <- split(seq_len(nrow(sub)), interaction(sub$system, sub$replicate, drop = TRUE))
  for (idx in groups) {
    seg <- sub[idx, , drop = FALSE]
    seg <- seg[order(seg$local_frame_index, seg$som_frame), , drop = FALSE]
    if (nrow(seg) < 2) {
      next
    }
    som_seg <- SOM
    som_seg$unit.classif <- SOM$unit.classif[seg$som_frame]
    if (!is.null(SOM$distances) && length(SOM$distances) >= max(seg$som_frame)) {
      som_seg$distances <- SOM$distances[seg$som_frame]
    }
    out <- out + comp.trans.mat(som_seg, start = 1)
  }
  out
}

make_population_som <- function(SOM, frame_map, systems = NULL) {
  sub <- frame_map
  if (!is.null(systems)) {
    sub <- sub[sub$system %in% systems, , drop = FALSE]
  }
  som_sub <- SOM
  som_sub$unit.classif <- SOM$unit.classif[sub$som_frame]
  if (!is.null(SOM$distances) && nrow(sub) > 0 && length(SOM$distances) >= max(sub$som_frame)) {
    som_sub$distances <- SOM$distances[sub$som_frame]
  }
  som_sub
}

write_transition_table <- function(transitions, label, out_tables, profile) {
  rows <- list()
  k <- 1
  for (i in seq_len(nrow(transitions))) {
    for (j in seq_len(ncol(transitions))) {
      value <- transitions[i, j]
      if (value > 0) {
        rows[[k]] <- data.frame(
          profile = profile,
          system = label,
          from_neuron = i,
          to_neuron = j,
          transitions = as.integer(value),
          edge_type = ifelse(i == j, "self", "between")
        )
        k <- k + 1
      }
    }
  }
  if (length(rows) == 0) {
    table <- data.frame(
      profile = character(),
      system = character(),
      from_neuron = integer(),
      to_neuron = integer(),
      transitions = integer(),
      edge_type = character()
    )
  } else {
    table <- do.call(rbind, rows)
  }
  outfile <- file.path(out_tables, paste0(profile, "_official_transition_edges_", label, ".csv"))
  write.csv(table, outfile, row.names = FALSE)
  invisible(outfile)
}

plot_network <- function(net, SOM, outfile, title, layout_type = c("som", "fr"), seed = 20260603) {
  layout_type <- match.arg(layout_type)
  if (layout_type == "som") {
    coords <- SOM$grid$pts
  } else {
    set.seed(seed)
    coords <- layout_with_fr(net, weights = NA)
  }
  png(outfile, width = 1800, height = 1400, res = 220)
  par(mar = c(1, 1, 4, 1))
  plot(
    net,
    edge.arrow.size = E(net)$width / 5,
    edge.curved = 0.17,
    edge.color = "black",
    vertex.label = "",
    layout = coords,
    main = title
  )
  dev.off()
}

main <- function() {
  args <- parse_args()
  dir.create(args$out_figures, recursive = TRUE, showWarnings = FALSE)
  dir.create(args$out_tables, recursive = TRUE, showWarnings = FALSE)

  SOM <- readRDS(file.path(args$run_dir, "SOM.rds"))
  frame_map <- read.csv(file.path(args$run_dir, "frame_mapping.csv"))
  neuron_clusters <- read.csv(file.path(args$run_dir, "neuron_clusters.csv"))
  SOM.hc <- neuron_clusters$cluster

  labels <- c("all", SYSTEM_ORDER)
  outputs <- c()
  edge_tables <- c()
  for (label in labels) {
    systems <- if (label == "all") NULL else label
    trans <- make_transition_matrix(SOM, frame_map, systems)
    som_pop <- make_population_som(SOM, frame_map, systems)
    edge_tables <- c(edge_tables, write_transition_table(trans, label, args$out_tables, args$profile))
    net <- matrix2graph(trans, som_pop, SOM.hc, COL_SCALE, diag = FALSE)

    outfile_som <- file.path(args$out_figures, paste0(args$profile, "_official_transition_network_", label, "_som_layout.png"))
    plot_network(
      net,
      SOM,
      outfile_som,
      paste(args$profile, label, "official SOMMD transition network: SOM layout"),
      "som",
      args$seed
    )
    outputs <- c(outputs, outfile_som)

    outfile_fr <- file.path(args$out_figures, paste0(args$profile, "_official_transition_network_", label, "_fr_layout.png"))
    plot_network(
      net,
      SOM,
      outfile_fr,
      paste(args$profile, label, "official SOMMD transition network: FR layout"),
      "fr",
      args$seed
    )
    outputs <- c(outputs, outfile_fr)
  }
  writeLines(c(outputs, edge_tables))
}

main()
