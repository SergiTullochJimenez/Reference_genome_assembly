# =============================================================================
# Genome Circos Plots — all 8 newly assembled Iberian snake genomes
# =============================================================================
# Six tracks per genome in 1-Mb windows (gene density, BUSCO assembly hits, repeat density, ONT coverage,
# Hi-C coverage, GC content), from the BED files written by 01.Circos_tracks.sh (R, circlize).
#
# Expected directory layout (BASE_DIR):
#
#   Circos_plot/
#     Coronella_austriaca/           <- one folder per species, containing:
#       chrom_sizes.txt
#       gc_1mb.bed
#       repeat_1mb.bed
#       gene_density_1mb.bed
#       busco_assembly_hits.bed
#       busco_annotation_hits.bed      (optional -- skipped if missing)
#       ont_cov_1mb.bed
#       hic_cov_1mb.bed
#     Coronella_girondica/
#     Hierophis_viridiflavus/
#     Macroprotodon_brevis/
#     Macroprotodon_mauritanicus/
#     Natrix_astreptophora/
#     Zamenis_longissimus/
#     Zamenis_scalaris/
#     Final_plots/                   <- created automatically if missing;
#                                        one circos_genome_<species>.pdf
#                                        per species lands here
#
# A species with missing/broken input files does NOT stop the other 7 --
# each species is wrapped in tryCatch, failures are logged and reported in
# a summary at the end.
# =============================================================================

BASE_DIR <- "/path/to/Circos_plot"                       # folder with one sub-folder of tracks per species
SEX_LISTS <- "/path/to/this/repo/analyses/macrosynteny"  # z_chromosome.txt / w_chromosome.txt
setwd(BASE_DIR)

suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
  library(RColorBrewer)
  library(scales)
})

# ── SPECIES LIST ──────────────────────────────────────────────────────────────
SPECIES_LIST <- c(
  "Coronella_austriaca",
  "Coronella_girondica",
  "Hierophis_viridiflavus",
  "Macroprotodon_brevis",
  "Macroprotodon_mauritanicus",
  "Natrix_astreptophora",
  "Zamenis_longissimus",
  "Zamenis_scalaris"
)

FINAL_PLOTS_DIR <- file.path(BASE_DIR, "Good")   # Z and W chromosomes drawn last
dir.create(FINAL_PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ── USER SETTINGS (same for every species) ────────────────────────────────────

MIN_SCAFFOLD     <- 5e6
ONT_TARGET_DEPTH <- NA
HIC_TARGET_DEPTH <- NA

# ── COLOUR PALETTE ────────────────────────────────────────────────────────────

PAL <- list(
  chrom_odd  = "#4A5240",
  chrom_even = "#8A9A7A",
  gc         = "#C8860A",
  repeat_col = "#7B3F00",
  gene       = "#2E6B8A",
  busco_c    = "#2D9B5A",
  busco_d    = "#F5C518",
  busco_f    = "#E07B39",
  busco_m    = "#C0392B",
  ont        = "#5B4A8A",
  hic        = "#1A7A6E"
)

busco_col <- c(
  Complete   = PAL$busco_c,
  Duplicated = PAL$busco_d,
  Fragmented = PAL$busco_f,
  Missing    = PAL$busco_m
)

# ── HELPERS ───────────────────────────────────────────────────────────────────

load_bed <- function(fname, data_dir, scaffolds, col_names, optional = FALSE) {
  f <- file.path(data_dir, fname)
  if (!file.exists(f)) {
    if (optional) { message("Optional file missing, skipping: ", fname); return(NULL) }
    stop("Missing input file: ", f)
  }
  dt <- fread(f, header = FALSE, col.names = col_names)
  dt <- dt[scaffold %in% scaffolds]
  if (nrow(dt) == 0) warning("No rows after scaffold filter: ", fname)
  dt
}

safe_ylim <- function(x, padding = 0.05) {
  mx <- max(x, na.rm = TRUE)
  if (!is.finite(mx) || mx == 0) mx <- 1
  c(0, mx * (1 + padding))
}

has_signal <- function(dt, col = "value") {
  !is.null(dt) && nrow(dt) > 0 && any(dt[[col]] > 0, na.rm = TRUE)
}

norm_cov  <- function(dt, target) {
  if (!is.na(target)) dt[, value := value / target]
  dt
}
winsorise <- function(x, p = 0.99) pmin(x, quantile(x, p, na.rm = TRUE))

# ── SCALE HELPER ──────────────────────────────────────────────────────────────
# Draws a proper y-axis on the left edge of the anchor sector's track cell.
# Call immediately after the circos.genomicTrack() it annotates.
#
# circos.yaxis() places the axis INSIDE the track, on the sector boundary,
# with tick marks pointing inward and labels sitting just outside -- exactly
# like a standard cartesian y-axis but bent to the circos geometry.
#
# Arguments:
#   track_index   : get.current.track.index() of the track just drawn
#   anchor_sector : sector to draw on (default: first/largest scaffold)
#   ylim          : c(min, max) matching the track's ylim
#   at            : explicit tick positions (overrides n_ticks if supplied)
#   n_ticks       : number of evenly-spaced ticks if `at` not supplied
#   labels        : explicit label strings (must match length of `at`)
#   label_suffix  : appended to auto-generated labels (e.g. "x", "%")
#   label_fmt     : sprintf format for auto-generated labels
#   col           : colour for axis line, ticks, and labels
#   cex           : label size

add_track_yaxis <- function(track_index,
                            anchor_sector,
                            ylim,
                            at           = NULL,
                            n_ticks      = 3,
                            labels       = NULL,
                            label_suffix = "",
                            label_fmt    = "%.1f",
                            col          = "grey30",
                            cex          = 0.40) {
  
  # Build tick positions
  if (is.null(at)) {
    at <- seq(ylim[1], ylim[2], length.out = n_ticks)
  }
  
  # Build tick labels
  if (is.null(labels)) {
    labels <- paste0(sprintf(label_fmt, at), label_suffix)
  }
  
  circos.yaxis(
    side         = "left",
    at           = at,
    labels       = labels,
    sector.index = anchor_sector,
    track.index  = track_index,
    col          = col,
    labels.col   = col,
    lwd          = 0.7,
    labels.cex   = cex,
    tick.length  = convert_x(0.5, "mm", sector.index = anchor_sector,
                             track.index = track_index)
  )
}

# =============================================================================
# PER-SPECIES PLOTTING FUNCTION
# =============================================================================

plot_circos_for_species <- function(species) {
  
  data_dir <- file.path(BASE_DIR, species)
  if (!dir.exists(data_dir)) {
    stop("No data folder found for ", species, " at ", data_dir)
  }
  out_file      <- file.path(FINAL_PLOTS_DIR,
                             paste0("circos_genome_", species, ".pdf"))
  species_label <- gsub("_", " ", species)
  
  # 1. LOAD & FILTER SCAFFOLDS ------------------------------------------------
  chrom_sizes_path <- file.path(data_dir, "chrom_sizes.txt")
  if (!file.exists(chrom_sizes_path)) stop("Missing input file: ", chrom_sizes_path)
  
  chrom_sizes <- fread(chrom_sizes_path,
                       col.names = c("scaffold", "length"), header = FALSE)
  chrom_sizes <- chrom_sizes[length >= MIN_SCAFFOLD][order(-length)]
  # Z (then W) go last: autosomes by size, then Z, then W (names from the macrosynteny z/w lists)
  sex_names <- function(f) {
    x <- fread(f, header = FALSE, sep = "\t", col.names = c("sp", "sc"), fill = TRUE)
    unlist(strsplit(x[sp == species]$sc, ","))
  }
  z_sc <- sex_names(file.path(SEX_LISTS, "z_chromosome.txt"))
  w_sc <- sex_names(file.path(SEX_LISTS, "w_chromosome.txt"))
  auto_sc <- setdiff(chrom_sizes$scaffold, c(z_sc, w_sc))
  new_order <- c(auto_sc, intersect(z_sc, chrom_sizes$scaffold), intersect(w_sc, chrom_sizes$scaffold))
  chrom_sizes <- chrom_sizes[match(new_order, scaffold)]
  cat(sprintf("[%s] last sectors: %s\n", species, paste(tail(new_order, 3), collapse = ", ")))
  scaffolds   <- chrom_sizes$scaffold
  anchor      <- scaffolds[1]   # largest scaffold — used for all scale axes
  
  cat(sprintf("[%s] Plotting %d scaffolds >= %.0f Mb  (total = %.2f Gb)\n",
              species, nrow(chrom_sizes), MIN_SCAFFOLD / 1e6,
              sum(chrom_sizes$length) / 1e9))
  
  # 2. LOAD TRACK DATA --------------------------------------------------------
  gc        <- load_bed("gc_1mb.bed",               data_dir, scaffolds, c("scaffold","start","end","value"))
  rep_den   <- load_bed("repeat_1mb.bed",           data_dir, scaffolds, c("scaffold","start","end","value"))
  genes     <- load_bed("gene_density_1mb.bed",     data_dir, scaffolds, c("scaffold","start","end","value"))
  busco_asm <- load_bed("busco_assembly_hits.bed",  data_dir, scaffolds, c("scaffold","start","end","status"))
  busco_ann <- load_bed("busco_annotation_hits.bed",data_dir, scaffolds, c("scaffold","start","end","status"),
                        optional = TRUE)
  ont       <- load_bed("ont_cov_1mb.bed",          data_dir, scaffolds, c("scaffold","start","end","value"))
  hic       <- load_bed("hic_cov_1mb.bed",          data_dir, scaffolds, c("scaffold","start","end","value"))
  
  # Normalise + winsorise coverage
  ont <- norm_cov(ont, ONT_TARGET_DEPTH); ont[, value := winsorise(value)]
  hic <- norm_cov(hic, HIC_TARGET_DEPTH); hic[, value := winsorise(value)]
  
  # Gene density: normalise to max; guard against all-zero
  max_genes <- max(genes$value, na.rm = TRUE)
  if (is.na(max_genes) || max_genes == 0) {
    warning("[", species, "] Gene density is all zero -- check GFF has 'gene' features on the plotted scaffolds")
    genes[, value_norm := 0]
  } else {
    genes[, value_norm := value / max_genes]
  }
  
  # Pre-compute ylims for tracks that need scales
  ont_ylim  <- safe_ylim(ont$value)
  hic_ylim  <- safe_ylim(hic$value)
  rep_ylim  <- c(0, 1)
  gene_ylim <- c(0, 1)
  
  # 3. BUILD CIRCOS LAYOUT ----------------------------------------------------
  pdf(out_file, width = 14, height = 14)
  on.exit(if (!is.null(dev.list())) dev.off(), add = TRUE)
  
  circos.par(
    gap.after               = c(rep(1.5, length(scaffolds) - 1), 6),
    start.degree            = 90,
    track.margin            = c(0.005, 0.01),
    cell.padding            = c(0, 0, 0, 0),
    points.overflow.warning = FALSE
  )
  
  circos.genomicInitialize(
    data         = chrom_sizes[, .(scaffold, start = 0L, end = length)],
    sector.names = scaffolds,
    plotType     = "axis"
  )
  
  for (sc in scaffolds) {
    circos.text(
      x            = mean(get.cell.meta.data("xlim", sector.index = sc)),
      y            = 1.8,
      labels       = gsub("_", " ", sc),
      sector.index = sc, track.index = 1,
      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
      cex = 0.42, col = "grey20"
    )
  }
  
  # TRACK 1 -- Gene density (normalised 0–1; axis labels show actual gene counts)
  circos.genomicTrack(
    genes[, .(scaffold, start, end, value = value_norm)],
    ylim = gene_ylim, track.height = 0.09, bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
                         ytop = value[[1]], ybottom = 0,
                         col = alpha(PAL$gene, 0.80), border = NA, ...)
    }
  )
  add_track_yaxis(
    track_index   = get.current.track.index(),
    anchor_sector = anchor,
    ylim          = gene_ylim,
    at            = c(0, 0.5, 1),
    labels        = c("0", sprintf("%.0f", max_genes * 0.5),
                      sprintf("%.0f", max_genes)),
    col           = PAL$gene
  )
  
  # TRACK 2 -- BUSCO assembly hits (dots; no continuous scale needed)
  circos.genomicTrack(
    busco_asm, ylim = c(0, 1), track.height = 0.04,
    bg.border = NA, bg.col = alpha("grey95", 0.3),
    panel.fun = function(region, value, ...) {
      col_vec <- busco_col[value[[1]]]
      y_vals  <- data.frame(y = rep(0.5, nrow(region)))
      circos.genomicPoints(region, y_vals,
                           col = col_vec, pch = 16, cex = 0.35, ...)
    }
  )
  
  # TRACK 3 -- BUSCO annotation hits (optional)
  has_ann <- !is.null(busco_ann) && nrow(busco_ann) > 0
  if (has_ann) {
    circos.genomicTrack(
      busco_ann, ylim = c(0, 1), track.height = 0.04,
      bg.border = NA, bg.col = alpha("grey95", 0.3),
      panel.fun = function(region, value, ...) {
        col_vec <- busco_col[value[[1]]]
        y_vals  <- data.frame(y = rep(0.5, nrow(region)))
        circos.genomicPoints(region, y_vals,
                             col = col_vec, pch = 16, cex = 0.35, ...)
      }
    )
  } else {
    message("[", species, "] BUSCO annotation: no busco_annotation_hits.bed (or 0 rows) -- track skipped")
  }
  
  # TRACK 4 -- Repeat density (0–100%; scale shows 0 / 50% / 100%)
  circos.genomicTrack(
    rep_den, ylim = rep_ylim, track.height = 0.09, bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
                         ytop = value[[1]], ybottom = 0,
                         col = alpha(PAL$repeat_col, 0.75), border = NA, ...)
    }
  )
  add_track_yaxis(
    track_index   = get.current.track.index(),
    anchor_sector = anchor,
    ylim          = rep_ylim,
    at            = c(0, 0.5, 1),
    labels        = c("0%", "50%", "100%"),
    col           = PAL$repeat_col
  )
  
  # TRACK 5 -- ONT coverage (scale shows 0 / mid / max depth)
  if (has_signal(ont)) {
    circos.genomicTrack(
      ont, ylim = ont_ylim, track.height = 0.09, bg.border = NA,
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value,
                            col = PAL$ont, lwd = 0.7, area = TRUE, border = NA, ...)
      }
    )
    add_track_yaxis(
      track_index   = get.current.track.index(),
      anchor_sector = anchor,
      ylim          = ont_ylim,
      at            = c(0, ont_ylim[2] / 2, ont_ylim[2]),
      labels        = c("0x", sprintf("%.0fx", ont_ylim[2] / 2),
                        sprintf("%.0fx", ont_ylim[2])),
      col           = PAL$ont
    )
  } else {
    message("[", species, "] ONT coverage: no signal detected -- track skipped")
  }
  
  # TRACK 6 -- HiC coverage (scale shows 0 / mid / max depth)
  if (has_signal(hic)) {
    circos.genomicTrack(
      hic, ylim = hic_ylim, track.height = 0.09, bg.border = NA,
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value,
                            col = PAL$hic, lwd = 0.7, area = TRUE, border = NA, ...)
      }
    )
    add_track_yaxis(
      track_index   = get.current.track.index(),
      anchor_sector = anchor,
      ylim          = hic_ylim,
      at            = c(0, hic_ylim[2] / 2, hic_ylim[2]),
      labels        = c("0x", sprintf("%.0fx", hic_ylim[2] / 2),
                        sprintf("%.0fx", hic_ylim[2])),
      col           = PAL$hic
    )
  } else {
    message("[", species, "] HiC coverage: no signal detected -- track skipped")
  }
  
  # TRACK 7 -- GC content (fixed 20–70% range; reference line at 40%)
  circos.genomicTrack(
    gc, ylim = c(0.2, 0.7), track.height = 0.09, bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value,
                          col = PAL$gc, lwd = 0.6, area = TRUE, border = NA,
                          baseline = 0.4, ...)
      circos.lines(
        x = c(0, get.cell.meta.data("xlim")[2]), y = c(0.4, 0.4),
        col = alpha("grey50", 0.5), lty = 2, lwd = 0.5)
    }
  )
  add_track_yaxis(
    track_index   = get.current.track.index(),
    anchor_sector = anchor,
    ylim          = c(0.2, 0.7),
    at            = c(0.2, 0.4, 0.7),
    labels        = c("20%", "40%", "70%"),
    col           = PAL$gc
  )
  
  # 4. LEGEND -----------------------------------------------------------------
  lgd_tracks <- list(
    list(label = "Gene density",   col = PAL$gene,       type = "r"),
    list(label = "BUSCO assembly", col = "grey40",       type = "p", pch = 16)
  )
  if (has_ann) lgd_tracks <- c(lgd_tracks,
                               list(list(label = "BUSCO annotation", col = "grey40", type = "p", pch = 16)))
  lgd_tracks <- c(lgd_tracks,
                  list(list(label = "Repeat density", col = PAL$repeat_col, type = "r")))
  if (has_signal(ont)) lgd_tracks <- c(lgd_tracks,
                                       list(list(label = "ONT coverage", col = PAL$ont, type = "l")))
  if (has_signal(hic)) lgd_tracks <- c(lgd_tracks,
                                       list(list(label = "HiC coverage", col = PAL$hic, type = "l")))
  lgd_tracks <- c(lgd_tracks,
                  list(list(label = "GC content", col = PAL$gc, type = "l")))
  
  legend("bottomleft",
         legend = sapply(lgd_tracks, `[[`, "label"),
         col    = sapply(lgd_tracks, `[[`, "col"),
         lwd    = sapply(lgd_tracks, function(x) if (x$type %in% c("l","r")) 4 else NA),
         pch    = sapply(lgd_tracks, function(x) if (x$type == "p") x$pch else NA),
         pt.cex = 1.2, bty = "n", cex = 0.65,
         title  = expression(bold("Tracks (outer -> inner)")), title.adj = 0
  )
  
  legend("bottomright",
         legend = names(busco_col), col = unlist(busco_col),
         pch = 16, pt.cex = 1.2, bty = "n", cex = 0.65,
         title = expression(bold("BUSCO status")), title.adj = 0
  )
  
  title(
    main     = paste0("Genome Assembly - ", species_label),
    sub      = sprintf("%d chromosomal scaffolds  |  1 Mb bins", nrow(chrom_sizes)),
    cex.main = 1.3, cex.sub = 0.7, col.sub = "grey50"
  )
  
  circos.clear()
  cat("[", species, "] Saved: ", out_file, "\n", sep = "")
}

# =============================================================================
# RUN ALL 8 SPECIES
# =============================================================================

results <- list()
for (species in SPECIES_LIST) {
  results[[species]] <- tryCatch({
    plot_circos_for_species(species)
    "OK"
  }, error = function(e) {
    message("[", species, "] FAILED: ", conditionMessage(e))
    paste("FAILED:", conditionMessage(e))
  })
}

cat("\n==================== SUMMARY ====================\n")
for (species in names(results)) {
  cat(sprintf("%-30s %s\n", species, results[[species]]))
}
n_ok <- sum(unlist(results) == "OK")
cat(sprintf("\n%d / %d species plotted successfully.\n", n_ok, length(SPECIES_LIST)))
cat("Plots saved to:", FINAL_PLOTS_DIR, "\n")

