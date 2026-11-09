#  ──────────────────────────────────────────────────────────────────────────
#  Time-scaled norovirus GII.4 phylogeny from a BEAST MCC tree.
#
#  Inputs
#    noro7.tree                     MCC tree (BEAST NEXUS, via read.beast)
#    iqtree_noro/noroiqtree.csv     metadata: taxa, capsid subtype, clinical_trait
#
#  ── 1. Load ───────────────────────────────────────────────────────────────
#
#    source("noro_phylo.R")
#    x <- phylo_load()                       # both files, default paths
#    print(x)
#
#  `x` holds the treedata, the metadata subset matching the tips, and the
#  most recent sampling date used to put the tree on a calendar axis.
#
#  ── 2. Plot ───────────────────────────────────────────────────────────────
#
#  One function draws every variant; the arguments switch layers on and off.
#
#    plot_phylo(x)                                  # tips by subtype
#    plot_phylo(x, heatmap = "clinical_trait")      # + trait heatmap strip
#    plot_phylo(x, node_colour = "location",        # + ancestral-state nodes
#               node_size = "location.prob")
#    plot_phylo(x, tip_colour = NULL,               # nodes only
#               node_colour = "location")
#
#  ── 3. Save ───────────────────────────────────────────────────────────────
#
#    p <- plot_phylo(x, heatmap = "clinical_trait",
#                    title = "Norovirus GII.4 VP1 phylogeny (n = 1183)")
#    ggsave("norotree.pdf", p, width = 18.4, height = 13.9, units = "in",
#           device = cairo_pdf)
#
#  ── 4. Different files ────────────────────────────────────────────────────
#
#    x <- phylo_load(tree = "~/noro8.tree",
#                    metadata = "~/iqtree_noro/noroiqtree_v2.csv",
#                    mrsd = "2025-03-01")
#
#  mrsd is the most recent sampling date and sets the calendar axis -- change
#  it whenever the sequence set changes, or every date on the tree is wrong.
#  Metadata needs a `taxa` column matching the tip labels; it is renamed to
#  `label` on load. Any other column can be used for colour or the heatmap.
#  Palettes live in the CONFIG block; a value with no entry is drawn grey.
#  ──────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(treeio)
library(ggtree)
library(ggnewscale)


# ── CONFIG ────────────────────────────────────────────────────────────────

PHYLO_MRSD   <- "2024-01-07"              # most recent sampling date
PHYLO_BREAKS <- seq(1970, 2025, 5)        # calendar axis ticks
PHYLO_NA     <- "#cccccc"                 # colour for unmatched values

SUBTYPE_PAL <- c(
  'Apeldoorn_2007'        = '#e6194b', 'Asia_2003'          = '#3cb44b',
  'Bristol_1993'          = '#f5a623', 'Camberwell_1994'    = '#4363d8',
  'Could not assign'      = '#aaaaaa', 'Den_Haag_2006b'     = '#f58231',
  'Farmington_Hills_2002' = '#911eb4', 'GII.3'              = 'black',
  'Hong_Kong_2019'        = '#f032e6', 'Hunter_2004'        = '#8bc34a',
  'Kaiso_2003'            = '#ff8a80', 'New_Orleans_2009'   = '#00897b',
  'Osaka_2007'            = '#ce93d8', 'Sydney_2012'        = '#795548',
  'Sydney_2016'           = '#b71c1c', 'US95_96'            = '#69f0ae',
  'Yerseke_2006a'         = '#3949ab'
)

TRAIT_PAL <- c(
  'acute'      = '#2166ac',
  'chronic'    = '#d6604d',
  'persistent' = '#1a9850',
  'unknown'    = 'grey30'
)

#  Which palette a variable gets. Add a row for any new annotation column.
PHYLO_PAL <- list(
  capsid.subtype = SUBTYPE_PAL,
  clinical_trait = TRAIT_PAL,
  location       = TRAIT_PAL
)

#  Legend titles, keyed the same way. Falls back to the column name.
PHYLO_LAB <- c(
  capsid.subtype = "Subtype / variant",
  clinical_trait = "Clinical trait",
  location       = "Ancestral state",
  location.prob  = "State probability",
  posterior      = "Posterior support"
)

#  Both lookups fall back rather than erroring, so a column with no entry
#  still plots -- with generated colours and the column name as its title.
phylo_pal <- function(v) if (v %in% names(PHYLO_PAL)) PHYLO_PAL[[v]] else NULL
phylo_lab <- function(v) if (v %in% names(PHYLO_LAB)) unname(PHYLO_LAB[v]) else v

scale_colour_for <- function(v, ...) {
  pal <- phylo_pal(v)
  if (is.null(pal))
    scale_colour_viridis_d(name = phylo_lab(v), na.value = PHYLO_NA, ...)
  else
    scale_colour_manual(name = phylo_lab(v), values = pal,
                        na.value = PHYLO_NA, ...)
}

scale_fill_for <- function(v, ...) {
  pal <- phylo_pal(v)
  if (is.null(pal))
    scale_fill_viridis_d(name = phylo_lab(v), na.value = PHYLO_NA, ...)
  else
    scale_fill_manual(name = phylo_lab(v), values = pal,
                      na.value = PHYLO_NA, ...)
}


# ── LOADING ───────────────────────────────────────────────────────────────

#' Read the MCC tree and its metadata.
#'
#' @param mrsd Most recent sampling date, "YYYY-MM-DD". Sets the calendar axis.
phylo_load <- function(tree     = "~/noro7.tree",
                       metadata = "~/iqtree_noro/noroiqtree.csv",
                       mrsd     = PHYLO_MRSD) {

  td <- read.beast(tree)
  tips <- td %>% as_tibble() %>% filter(!is.na(label)) %>% pull(label)

  meta <- read.csv(metadata, stringsAsFactors = FALSE) %>%
    rename(label = taxa) %>%
    mutate(across(everything(), trimws)) %>%
    filter(label %in% tips)

  #  gheatmap() matches its data to tips by row name, not by column.
  rownames(meta) <- meta$label

  missing <- setdiff(tips, meta$label)
  if (length(missing))
    warning(length(missing), " tip(s) have no metadata row, e.g. ",
            paste(head(missing, 3), collapse = ", "), call. = FALSE)

  structure(list(tree = td, meta = meta, mrsd = mrsd,
                 mrsd_dec = lubridate::decimal_date(as.Date(mrsd))),
            class = "norophylo")
}

print.norophylo <- function(x, ...) {
  cat("<norophylo>\n")
  cat("  tips    :", ape::Ntip(as.phylo(x$tree)), "\n")
  cat("  metadata:", nrow(x$meta), "rows;",
      paste(setdiff(names(x$meta), "label"), collapse = ", "), "\n")
  cat("  node data:", paste(setdiff(names(x$tree@data), "node"),
                            collapse = ", "), "\n")
  cat("  mrsd    :", x$mrsd, sprintf("(%.4f)", x$mrsd_dec), "\n")
  invisible(x)
}

#' Node table with calendar dates, handy for picking clades to highlight.
phylo_nodes <- function(x) {
  x$tree %>% as_tibble() %>%
    transmute(node, label,
              date = x$mrsd_dec - as.numeric(height),
              across(any_of(c("location", "location.prob", "posterior"))))
}


# ── THEME ─────────────────────────────────────────────────────────────────

#  Two gotchas, both of which silently remove the time axis:
#   * theme_tree2() takes no base_size -- it forwards unknown arguments to
#     theme(), where `base_size` is not a valid element. Set it via `text`.
#   * returning `theme_tree2() + theme(...)` drops the theme's "complete"
#     flag, so ggtree's own axis-blanking theme wins and axis.text.x comes
#     out element_blank. Returning a LIST applies them in order instead.
theme_phylo <- function(base_size = 11) {
  list(
    theme_tree2(),
    theme(
      text            = element_text(size = base_size),
      axis.text.x     = element_text(size = rel(0.9), colour = "grey20"),
      axis.line.x     = element_line(colour = "grey20", linewidth = 0.3),
      axis.ticks.x    = element_line(colour = "grey20", linewidth = 0.3),
      legend.position = "right",
      legend.box      = "vertical",
      legend.key      = element_blank(),
      legend.text     = element_text(size = rel(1)),
      legend.title    = element_text(size = rel(1.1), face = "bold"),
      plot.title      = element_text(size = rel(1.4), face = "bold"),
      plot.caption    = element_text(size = rel(0.9), colour = "grey30")
    )
  )
}


# ── THE PLOT ──────────────────────────────────────────────────────────────

#' Time-scaled tree with optional tip colours, ancestral-state nodes and a
#' metadata heatmap strip.
#'
#' @param tip_colour  Metadata column for tip points, or NULL for none.
#' @param node_colour Node-data column for internal-node points, or NULL.
#' @param node_size   Node-data column scaling those points (e.g. a posterior).
#' @param heatmap     Metadata column drawn as a strip beside the tree, or NULL.
#' @param node_dates  TRUE writes the calendar date beside each internal node.
plot_phylo <- function(x,
                       tip_colour   = "capsid.subtype",
                       node_colour  = NULL,
                       node_size    = NULL,
                       heatmap      = NULL,
                       node_dates   = FALSE,
                       tip_size     = 0.75,
                       node_range   = c(0.5, 3),
                       heatmap_width = 0.05,
                       heatmap_offset = 0.01,
                       breaks       = PHYLO_BREAKS,
                       base_size    = 11,
                       title        = NULL) {

  p <- ggtree(x$tree, ladderize = TRUE, mrsd = x$mrsd) %<+% x$meta

  #  Internal nodes first, so tips are drawn on top of them.
  if (!is.null(node_colour)) {
    node_aes <- if (is.null(node_size)) {
      aes(colour = .data[[node_colour]], subset = !isTip)
    } else {
      aes(colour = .data[[node_colour]], subset = !isTip,
          size = as.numeric(.data[[node_size]]))
    }
    p <- p + geom_nodepoint(node_aes) +
      scale_colour_for(node_colour) +
      guides(colour = guide_legend(order = 1,
                                   override.aes = list(size = 4)))
    if (!is.null(node_size))
      p <- p + scale_size_continuous(name = phylo_lab(node_size),
                                     range = node_range)
    if (!is.null(tip_colour)) p <- p + new_scale_colour()
  }

  if (!is.null(tip_colour)) {
    p <- p + geom_tippoint(aes(colour = .data[[tip_colour]]), size = tip_size) +
      scale_colour_for(tip_colour, guide = guide_legend(ncol = 1)) +
      guides(colour = guide_legend(order = 2,
                                   override.aes = list(size = 4)))
  }

  if (node_dates)
    p <- p + geom_nodelab(aes(label = round(x$mrsd_dec - as.numeric(height), 1),
                              subset = !isTip),
                          size = 1, colour = "black", hjust = -0.3)

  #  gheatmap brings its own fill scale, so give it a fresh one.
  if (!is.null(heatmap)) {
    p <- p + new_scale_fill()
    p <- gheatmap(p, x$meta[, heatmap, drop = FALSE],
                  offset = heatmap_offset, width = heatmap_width,
                  colnames = FALSE, color = NA) +
      scale_fill_for(heatmap)
  }

  p +
    scale_x_continuous(breaks = breaks) +
    labs(caption = "Year", title = title) +
    theme_phylo(base_size)
}
