library(treeio)
library(ggtree)
library(ggnewscale)
library(ggplot2)
library(dplyr)
library(ape)
library(phangorn)
library(tibble)

#' Plot an annotated norovirus GII.4 VP1 tree (ML)
#'
#' @param treefile   path to Newick tree file
#' @param metadata   path to metadata CSV, or a data.frame already in memory
#' @param outgroup   tip label to root on and drop (NULL = just midpoint root, keep all tips)
#' @param label_col  metadata column that matches tree tip labels (default "taxa")
#' @param subtype_col metadata column for subtype/variant tip color (default "capsid.subtype")
#' @param trait_col  metadata column for the clinical trait heatmap strip (default "clinical_trait")
#' @param tree_title plot title; auto-generated with tip count if NULL
#' @param out_prefix output filename prefix (no extension)
#' @param save       write PDF (+ PNG) to disk
plot_noro_tree <- function(treefile,
                           metadata,
                           outgroup       = "U02030.1_851-2497",
                           label_col      = "taxa",
                           subtype_col    = "capsid.subtype",
                           trait_col      = "clinical_trait",
                           tree_title     = NULL,
                           heatmap_offset = 0.01,
                           heatmap_width  = 0.05,
                           out_prefix     = "norotree",
                           width          = 18.38,
                           height         = 13.94,
                           save           = TRUE) {
  
  # ── palettes ────────────────────────────────────────────────────────────
  subtype_pal <- c(
    'Apeldoorn_2007'        = '#e6194b', 'Asia_2003'             = '#3cb44b',
    'Bristol_1993'          = '#f5a623', 'Camberwell_1994'       = '#4363d8',
    'Could not assign'      = '#aaaaaa', 'Den_Haag_2006b'        = '#f58231',
    'Farmington_Hills_2002' = '#911eb4', 'GII.3'                 = 'black',
    'Hong_Kong_2019'        = '#f032e6', 'Hunter_2004'           = '#8bc34a',
    'Kaiso_2003'            = '#ff8a80', 'New_Orleans_2009'      = '#00897b',
    'Osaka_2007'            = '#ce93d8', 'Sydney_2012'           = '#795548',
    'Sydney_2016'           = '#b71c1c', 'US95_96'               = '#69f0ae',
    'Yerseke_2006a'         = '#3949ab'
  )
  trait_pal <- c(
    'acute'      = '#0072B2',
    'chronic'    = '#D55E00',
    'persistent' = '#009E73',
    'unknown'    = '#7F7F7F'
  )
  
  # ── tree ────────────────────────────────────────────────────────────────
  tr <- read.newick(treefile)
  tr <- midpoint(tr)
  if (!is.null(outgroup)) {
    tr <- root(tr, outgroup = outgroup, edgelabel = TRUE)
    tr <- drop.tip(tr, outgroup)
  }
  n <- length(tr$tip.label)
  
  # ── metadata ────────────────────────────────────────────────────────────
  meta <- if (is.character(metadata)) read.csv(metadata, stringsAsFactors = FALSE) else metadata
  meta <- meta %>%
    rename(label = !!label_col) %>%
    mutate(across(everything(), trimws))
  
  if (is.null(tree_title)) {
    tree_title <- sprintf(
      "Norovirus GII.4 VP1 phylogeny with clinical traits and subtypes (n = %d)", n
    )
  }
  
  # ── base tree + subtype tip points ─────────────────────────────────────
  p <- ggtree(tr) %<+% meta +
    geom_tippoint(aes(color = .data[[subtype_col]]), size = 0.5) +
    scale_color_manual(name = "Subtype/Variant", values = subtype_pal,
                       na.value = "#cccccc", guide = guide_legend(ncol = 1)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_tree() +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 18),
      legend.title    = element_text(size = 20, face = "bold")
    )
  
  # ── clinical trait heatmap strip ───────────────────────────────────────
  trait_df <- meta %>%
    select(label, all_of(trait_col)) %>%
    column_to_rownames("label")
  
  final_plot <- gheatmap(
    p + new_scale_fill(), trait_df,
    offset = heatmap_offset, width = heatmap_width,
    colnames = FALSE, color = NA
  ) +
    scale_fill_manual(name = "Clinical Trait", values = trait_pal, na.value = "black") +
    theme_tree() +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 18),
      legend.title    = element_text(size = 20, face = "bold"),
      plot.title      = element_text(size = 24, face = "bold")
    ) +
    ggtitle(tree_title)
  
  if (save) {
    ggsave(paste0(out_prefix, ".pdf"), final_plot, width = width, height = height,
           units = "in", device = cairo_pdf)
    ggsave(paste0(out_prefix, ".png"), final_plot, width = width, height = height,
           units = "in", dpi = 300)
  }
  
  final_plot
}

# ── usage ─────────────────────────────────────────────────────────────────
p <- plot_noro_tree(
  treefile = "~/iqtree_noro/noro2_aligned.fasta.treefile",
  metadata = "~/iqtree_noro/noroiqtree.csv"
)
p


