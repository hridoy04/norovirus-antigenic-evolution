#  ──────────────────────────────────────────────────────────────────────────
#  Norovirus GII.4 VP1 mutation figures.
#
#  Two input files, both in the home directory by default:
#    noro4_acpu.csv     position, mut_no, group
#    noro_aa_probs.csv  clinical_trait, pos, alt, probability
#
#  Usage
#  ─────
#    source("noro_mutation.R")
#    x <- noro_load(n_seq = c(acute = 605, chronic = 353,
#                             persistent = 9, unknown = 37))
#    noro_states(x)                                  # what is available
#
#    plot_site_freq(x, states = c("acute", "chronic"))
#    plot_aa_probs (x, states = c("acute", "chronic"))
#    plot_aa_legend()
#
#  `states` everywhere accepts:
#    NULL                                  every state in the file
#    c("chronic", "acute")                 a subset, drawn in this order
#    c(Acute = "acute", Chronic = "chronic")   subset, relabelled
#
#  Everything meant for tweaking is in the CONFIG block.
#  ──────────────────────────────────────────────────────────────────────────

library(tidyverse)


# ── CONFIG ────────────────────────────────────────────────────────────────

#  Antigenic sites, named by letter. Edit these and every figure follows.
NORO_SITES <- list(
  A = c(294:298, 368, 372, 373),
  B = c(333, 389),
  C = c(339:341, 375:378),
  D = c(393:398),
  E = c(407:414),
  G = c(352, 355:357, 359, 364)
)

#  One hue per site, taken in order from a documented CVD-checked categorical
#  order (worst adjacent-pair CVD dE 9.1, normal-vision dE 19.6). Positions
#  outside every site fall back to NORO_OTHER.
NORO_SITE_PAL <- c(
  A = "#2a78d6",   # blue
  B = "#eb6834",   # orange
  C = "#1baf7a",   # aqua
  D = "#eda100",   # yellow
  E = "#e87ba4",   # magenta
  G = "#6250d6"    # violet
)
NORO_OTHER <- "#c9c8c3"

#  How a site is written on the figure.
noro_site_label <- function(s) if_else(s == "other", "Outside sites",
                                       paste("Site", s))

NORO_INK <- c(primary = "#0b0b0b", secondary = "#3d3c3a",
              muted = "#7a7974", rule = "#2b2a28", grid = "#e6e5e1")

#  Typeface for the figures. "" uses the device default; set to a font you
#  have installed (e.g. "Helvetica", "Arial", "Source Sans Pro") for print.
NORO_FONT <- ""

#  VP1 regions available to the `region` argument.
NORO_REGIONS <- list(
  p2  = c(279, 406),
  p1  = c(226, 278),
  s   = c(1, 225),
  all = c(1, 540)
)

#  Amino-acid colours, grouped by chemical property.
AA_PALETTE <- c(
  # Hydrophobic (blues)
  A = "#C8D7FF", V = "#A0C0FF", I = "#78A8FF", L = "#5090FF",
  M = "#2878FF", F = "#0060FF", W = "#0048D1",
  # Polar (greens)
  S = "#C8FFC8", T = "#96FF96", N = "#64FF64", Q = "#32CD32",
  # Positive (reds)
  R = "#FFB4B4", K = "#FF6E6E", H = "#FF0000",
  # Negative (purples)
  D = "#E6BEFF", E = "#9B30FF",
  # Special
  C = "#F08080", G = "#F09048", P = "#C0C000", Y = "#15A4A4",
  # Gaps / unknown
  `-` = "#E0E0E0", X = "#B0B0B0", `?` = "#D0D0D0", `*` = "#A0A0A0",
  B = "#C0C0C0", Z = "#C0C0C0", J = "#C0C0C0"
)


# ── LOADING ───────────────────────────────────────────────────────────────

#  Map a position to its antigenic site, or NA.
noro_site_of <- function(pos, sites = NORO_SITES) {
  out <- rep(NA_character_, length(pos))
  for (s in names(sites)) out[pos %in% sites[[s]]] <- s
  out
}

#' Read both CSVs and return a `noro` object.
#'
#' @param n_seq Named vector of sequences per group, e.g.
#'   c(acute = 605, chronic = 353, persistent = 9, unknown = 37).
#'   This is the denominator for `freq`. If NULL it falls back to the largest
#'   `mut_no` seen in each group, which is only a lower bound -- pass the real
#'   counts before publishing anything.
noro_load <- function(acpu  = "~/noro4_acpu.csv",
                      probs = "~/noro_aa_probs.csv",
                      n_seq = NULL) {

  acpu <- read_csv(acpu, show_col_types = FALSE) %>%
    mutate(pos = as.integer(as.character(position)),
           site = noro_site_of(pos))

  probs <- read_csv(probs, show_col_types = FALSE) %>%
    select(-any_of("...1")) %>%
    mutate(pos = as.integer(pos),
           site = noro_site_of(pos))

  if (is.null(n_seq)) {
    n_seq <- acpu %>% group_by(group) %>%
      summarise(n = max(mut_no), .groups = "drop") %>% deframe()
    message("n_seq not supplied; using max(mut_no) per group as the denominator: ",
            paste(names(n_seq), n_seq, sep = "=", collapse = ", "),
            "\n  These are lower bounds -- pass the real sequence counts for ",
            "publication figures.")
  }

  acpu <- acpu %>%
    mutate(total_no = unname(n_seq[group]),
           freq     = mut_no / total_no)

  structure(list(acpu = acpu, probs = probs, n_seq = n_seq),
            class = "noro")
}

print.noro <- function(x, ...) {
  cat("<noro>\n")
  cat("  acpu :", nrow(x$acpu), "rows,", n_distinct(x$acpu$pos), "positions\n")
  cat("  probs:", nrow(x$probs), "rows,", n_distinct(x$probs$pos), "positions\n")
  cat("  states:", paste(noro_states(x), collapse = ", "), "\n")
  cat("  n_seq :", paste(names(x$n_seq), x$n_seq, sep = "=", collapse = ", "), "\n")
  invisible(x)
}

#' The clinical states present in both tables.
noro_states <- function(x) {
  intersect(sort(unique(x$acpu$group)), sort(unique(x$probs$clinical_trait)))
}

#  Turn the `states` argument into c(level = label), validating as it goes.
#  Accepts NULL (all), a bare vector, or a named vector of labels.
noro_resolve_states <- function(states, available) {
  if (is.null(states)) {
    out <- setNames(str_to_sentence(available), available)
  } else {
    lv  <- unname(states)
    lab <- if (is.null(names(states))) str_to_sentence(lv) else names(states)
    bad <- setdiff(lv, available)
    if (length(bad))
      stop("Unknown clinical state(s): ", paste(bad, collapse = ", "),
           ".\n  Available: ", paste(available, collapse = ", "), call. = FALSE)
    out <- setNames(lab, lv)
  }
  out
}

#  Resolve `region` to c(min, max).
noro_region <- function(region) {
  if (is.numeric(region) && length(region) == 2) return(region)
  NORO_REGIONS[[match.arg(tolower(region), names(NORO_REGIONS))]]
}

#' The positions a figure will show -- shared by both plots, so the amino-acid
#' panel lines up with the frequency panel by default.
noro_positions <- function(x,
                           region     = "p2",
                           states     = NULL,
                           sites_only = FALSE,
                           min_freq   = 0) {
  r <- noro_region(region)
  d <- x$acpu %>% filter(pos >= r[1], pos <= r[2], freq >= min_freq)
  if (!is.null(states)) d <- d %>% filter(group %in% unname(states))
  if (sites_only)       d <- d %>% filter(!is.na(site))
  sort(unique(d$pos))
}


# ── THEME ─────────────────────────────────────────────────────────────────

#' Print-oriented theme: hairline axes, no grid, real black text, generous
#' panel spacing, legend as a single row beneath the plot.
theme_noro <- function(base_size = 9, base_family = NORO_FONT) {
  half <- base_size / 2
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text             = element_text(colour = NORO_INK[["primary"]]),
      plot.title       = element_text(face = "bold", size = rel(1.15),
                                      margin = margin(b = half * 0.6)),
      plot.subtitle    = element_text(colour = NORO_INK[["secondary"]],
                                      size = rel(0.95),
                                      margin = margin(b = half * 1.2)),
      plot.caption     = element_text(colour = NORO_INK[["muted"]],
                                      size = rel(0.8), hjust = 0,
                                      margin = margin(t = half)),
      axis.title.x     = element_text(size = rel(1), margin = margin(t = half)),
      axis.title.y     = element_text(size = rel(1), margin = margin(r = half)),
      axis.text        = element_text(colour = NORO_INK[["secondary"]],
                                      size = rel(0.85)),
      axis.line        = element_line(colour = NORO_INK[["rule"]],
                                      linewidth = 0.3),
      axis.ticks       = element_line(colour = NORO_INK[["rule"]],
                                      linewidth = 0.3),
      axis.ticks.length = unit(2.5, "pt"),
      panel.grid       = element_blank(),
      panel.spacing    = unit(half * 1.4, "pt"),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "horizontal",
      legend.margin    = margin(t = half * 0.5),
      legend.key.size  = unit(base_size * 0.9, "pt"),
      legend.title     = element_text(size = rel(0.95),
                                      colour = NORO_INK[["secondary"]]),
      legend.text      = element_text(size = rel(0.9),
                                      colour = NORO_INK[["secondary"]]),
      strip.background = element_blank(),
      strip.text.y.right = element_text(face = "bold", size = rel(1),
                                        angle = -90,
                                        colour = NORO_INK[["primary"]]),
      strip.text       = element_text(face = "bold", size = rel(1),
                                      colour = NORO_INK[["primary"]]),
      plot.margin      = margin(half, half, half, half),
      plot.title.position   = "plot",
      plot.caption.position = "plot"
    )
}


# ── PLOT 1: mutation frequency, bars coloured by antigenic site ───────────

#' Mutation frequency per position, one panel per clinical state.
#'
#' @param states See the header: NULL, a subset, or a named vector of labels.
#' @param region "p2" (default), "p1", "s", "all", or c(min, max).
#' @param sites_only TRUE keeps only antigenic-site positions and puts them on
#'   a discrete axis.
#' @param label_sites Off by default -- the bar colours and the legend already
#'   carry the site. TRUE writes the site name above the top panel as well.
plot_site_freq <- function(x,
                           states      = NULL,
                           region      = "p2",
                           sites_only  = FALSE,
                           label_sites = FALSE,
                           free_y      = TRUE,
                           pal         = NORO_SITE_PAL,
                           base_size   = 9,
                           title       = NULL,
                           subtitle    = NULL) {

  st <- noro_resolve_states(states, noro_states(x))
  r  <- noro_region(region)

  d <- x$acpu %>%
    filter(pos >= r[1], pos <= r[2], group %in% names(st))
  if (sites_only) d <- d %>% filter(!is.na(site))
  if (!nrow(d)) stop("No rows left after filtering.", call. = FALSE)

  d <- d %>%
    mutate(group  = factor(unname(st[group]), levels = unname(st)),
           site_f = factor(replace_na(site, "other"),
                           levels = c(names(pal), "other")))

  fills  <- c(pal, other = NORO_OTHER)
  shown  <- intersect(names(pal), unique(as.character(d$site_f)))

  p <- ggplot(d, aes(y = freq, fill = site_f))

  if (sites_only) {
    p <- p + geom_col(aes(x = factor(pos, levels = sort(unique(pos)))),
                      width = 0.72) +
      labs(x = "Antigenic-site position")
  } else {
    p <- p + geom_col(aes(x = pos), width = 0.9) +
      scale_x_continuous(breaks = scales::breaks_width(10),
                         expand = expansion(mult = 0.012)) +
      labs(x = "VP1 codon position")
  }

  #  One label per site, over the middle of its widest contiguous run, so a
  #  scattered site like G is named once rather than at every run.
  if (label_sites && !sites_only) {
    runs <- imap_dfr(NORO_SITES, ~ tibble(site = .y, pos = .x)) %>%
      filter(pos >= r[1], pos <= r[2]) %>%
      arrange(site, pos) %>%
      group_by(site) %>%
      mutate(run = cumsum(c(1, diff(pos) != 1))) %>%
      group_by(site, run) %>%
      summarise(x = mean(pos), len = n(), .groups = "drop_last") %>%
      slice_max(len, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(lab = noro_site_label(site),
             group = factor(levels(d$group)[1], levels = levels(d$group)))
    p <- p + geom_text(data = runs, aes(x = x, y = Inf, label = lab),
                       inherit.aes = FALSE, vjust = 1.3,
                       size = base_size * 0.30, fontface = "bold",
                       family = NORO_FONT, colour = NORO_INK[["secondary"]])
  }

  p +
    scale_fill_manual(name   = "Antigenic site",
                      values = fills,
                      breaks = shown,
                      labels = noro_site_label(shown),
                      na.translate = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, if (label_sites) 0.20 else 0.05)),
                       labels = scales::label_percent(accuracy = 1)) +
    facet_wrap(~ group, ncol = 1,
               scales = if (free_y) "free_y" else "fixed",
               strip.position = "right") +
    labs(
      y = "Sequences mutated at this position",
      title = title,
      subtitle = subtitle
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_noro(base_size) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                     size = rel(0.72)))
}


# ── PLOT 2: amino-acid composition at those same positions ────────────────

#' Amino-acid probability per position, stacked, one panel per clinical state.
#'
#' @param positions Default: the antigenic-site positions plot_site_freq()
#'   would draw for the same `region` and `states`, so the two figures share
#'   an x-axis.
#' @param drop_fixed TRUE hides positions where one residue is at probability
#'   1 in every state shown -- i.e. nothing varies there.
plot_aa_probs <- function(x,
                          positions  = NULL,
                          states     = NULL,
                          region     = "p2",
                          sites_only = TRUE,
                          drop_fixed = TRUE,
                          label_aa   = TRUE,
                          site_strip = TRUE,
                          base_size  = 9,
                          title      = NULL,
                          subtitle   = NULL) {

  st <- noro_resolve_states(states, noro_states(x))

  if (is.null(positions)) {
    positions <- noro_positions(x, region = region, states = names(st),
                                sites_only = sites_only)
  }

  d <- x$probs %>%
    filter(pos %in% positions, clinical_trait %in% names(st))
  if (!nrow(d)) stop("No rows left after filtering.", call. = FALSE)

  if (drop_fixed) {
    varying <- d %>% group_by(pos) %>%
      summarise(varies = any(probability < 1), .groups = "drop") %>%
      filter(varies) %>% pull(pos)
    d <- d %>% filter(pos %in% varying)
    if (!nrow(d))
      stop("Every selected position is invariant; set drop_fixed = FALSE.",
           call. = FALSE)
  }

  d <- d %>%
    mutate(clinical_trait = factor(unname(st[clinical_trait]),
                                   levels = unname(st)),
           pos_f = factor(pos, levels = sort(unique(pos))),
           site  = replace_na(site, "other"))

  p <- ggplot(d, aes(pos_f, probability, fill = alt)) +
    geom_col(width = 0.82, colour = "white", linewidth = 0.2)

  #  Residue letters inside segments tall enough to hold one.
  if (label_aa) {
    p <- p + geom_text(aes(label = if_else(probability >= 0.12, alt, "")),
                       position = position_stack(vjust = 0.5),
                       size = base_size * 0.26, family = NORO_FONT,
                       colour = NORO_INK[["primary"]])
  }

  #  The letters on the bars are the key; the residue legend is only needed
  #  when they are switched off (plot_aa_legend() is the standalone version).
  p <- p + scale_fill_manual(name = "Residue", values = AA_PALETTE,
                             guide = if (label_aa) "none" else "legend")

  #  A colour strip under the bars carrying the antigenic site, so this panel
  #  and plot_site_freq() can be read against each other.
  if (site_strip) {
    strip <- d %>% distinct(pos_f, site)
    shown <- intersect(c(names(NORO_SITE_PAL), "other"), unique(strip$site))
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_tile(data = strip, inherit.aes = FALSE,
                aes(x = pos_f, y = -0.055, fill = site),
                height = 0.045, width = 0.82) +
      scale_fill_manual(name   = "Antigenic site",
                        values = c(NORO_SITE_PAL, other = NORO_OTHER),
                        breaks = shown,
                        labels = noro_site_label(shown)) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  }

  p +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       labels = scales::label_percent(accuracy = 1)) +
    facet_wrap(~ clinical_trait, ncol = 1, strip.position = "right") +
    labs(x = "Antigenic-site position", y = "Amino-acid frequency",
         title = title, subtitle = subtitle) +
    theme_noro(base_size) +
    theme(legend.position = if (site_strip) "bottom" else "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}


# ── Standalone amino-acid colour key ──────────────────────────────────────

plot_aa_legend <- function(pal = AA_PALETTE, base_size = 9) {
  groups <- tribble(
    ~aa, ~group,
    "A","Hydrophobic","V","Hydrophobic","I","Hydrophobic","L","Hydrophobic",
    "M","Hydrophobic","F","Hydrophobic","W","Hydrophobic",
    "S","Polar","T","Polar","N","Polar","Q","Polar",
    "R","Positive","K","Positive","H","Positive",
    "D","Negative","E","Negative",
    "C","Special","G","Special","P","Special","Y","Special",
    "-","Gap/unknown","X","Gap/unknown","?","Gap/unknown","*","Gap/unknown"
  ) %>%
    group_by(group) %>% mutate(order = row_number()) %>% ungroup() %>%
    mutate(group = fct_rev(fct_inorder(group)))

  ggplot(groups, aes(order, group, fill = aa)) +
    geom_tile(colour = "white", width = 0.9, height = 0.9) +
    geom_text(aes(label = aa), size = base_size * 0.30, fontface = "bold",
              family = NORO_FONT) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_noro(base_size) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_text(face = "bold",
                                     colour = NORO_INK[["primary"]]))
}





# ── EXAMPLES ──────────────────────────────────────────────────────────────
#  Commented out so that source() only defines functions. Copy a block into
#  the console, or un-comment the ones you use.

## ---- Setup ---------------------------------------------------------------
#  You source THIS file from the console or your analysis script:
#      source("noro_mutation.R")
#  Never un-comment a source() of this file here -- it would source
#  itself forever and die with "C stack usage is too close to the limit".
# library(ggnewscale)                       # needed by the site strip
# 
# x <- noro_load(n_seq = c(acute = 605, chronic = 353,
#                          persistent = 9, unknown = 37))
# print(x)                                  # states + denominators
# noro_states(x)                            # -> "acute" "chronic" "persistent" "unknown"


## ---- Choosing clinical states --------------------------------------------
# All states in the file, one panel each:
# plot_site_freq(x)
#
# Just two, in the order you type them:
# plot_site_freq(x, states = c("acute", "chronic"))
# plot_site_freq(x, states = c("chronic", "acute"))     # panels swap
#
# The small groups. free_y = FALSE puts them on a shared scale, which is the
# honest view when the group sizes differ a lot:
# plot_site_freq(x, states = c("persistent", "unknown"), free_y = FALSE)
#
# One state on its own:
# plot_site_freq(x, states = "chronic")
#
# Rename the panels for the figure without touching the data. Names become
# the labels, values must match the file:
# plot_site_freq(x, states = c("Acute infection"   = "acute",
#                              "Chronic infection" = "chronic"))
#
# The same argument works identically on the other plot:
# plot_aa_probs(x, states = c("persistent", "unknown"))
#
# A typo tells you what is valid rather than drawing an empty panel:
# plot_site_freq(x, states = "acutte")
# # Error: Unknown clinical state(s): acutte.
# #   Available: acute, chronic, persistent, unknown


## ---- Updated CSVs --------------------------------------------------------
# Point noro_load() at the new files. Nothing else changes -- every plot
# takes x, so they all follow:
#
# x2 <- noro_load(acpu  = "~/noro4_acpu_v2.csv",
#                 probs = "~/noro_aa_probs_v2.csv",
#                 n_seq = c(acute = 710, chronic = 402,
#                           persistent = 12, unknown = 41))
# plot_site_freq(x2, states = c("acute", "chronic"))
#
# Don't know the new sequence counts yet? Omit n_seq for a quick look. It
# falls back to max(mut_no) per group and says so:
# x2 <- noro_load("~/noro4_acpu_v2.csv", "~/noro_aa_probs_v2.csv")
#
# Always check a new load before trusting the y-axis. Any max_freq above 1
# is impossible and means n_seq is too small for that group:
# print(x2)
# x2$acpu %>% count(group)
# x2$acpu %>% group_by(group) %>% summarise(max_freq = max(freq))
#
# If the update introduces a new state, it just appears -- no code change:
# noro_states(x2)
# plot_site_freq(x2, states = c("acute", "reinfection"))
# # ...but remember to add it to n_seq, or its bars come out NA and vanish.


## ---- Regions and positions -----------------------------------------------
# plot_site_freq(x, region = "all")                    # whole VP1
# plot_site_freq(x, region = c(330, 420))              # custom window
# plot_site_freq(x, states = c("acute","chronic"), sites_only = TRUE)
#
# By default plot_aa_probs() uses the positions plot_site_freq() would draw.
# Override with your own list:
# plot_aa_probs(x, positions = c(294, 297, 340, 368, 373, 377))
#
# Keep positions where a single residue is fixed at 100%:
# plot_aa_probs(x, drop_fixed = FALSE)


## ---- Saving --------------------------------------------------------------
# p1 <- plot_site_freq(x, states = c("acute", "chronic"))
# p2 <- plot_aa_probs (x, states = c("acute", "chronic"))
#
# ggsave("fig_sitefreq.pdf", p1, width = 180, height = 95,  units = "mm",
#        device = cairo_pdf)
# ggsave("fig_aaprobs.pdf",  p2, width = 180, height = 100, units = "mm",
#        device = cairo_pdf)
# ggsave("fig_aalegend.pdf", plot_aa_legend(), width = 110, height = 55,
#        units = "mm", device = cairo_pdf)
#
# Single-column version -- narrower, so bump the type up a little:
# ggsave("fig_sitefreq_1col.pdf",
#        plot_site_freq(x, states = c("acute","chronic"), base_size = 11),
#        width = 90, height = 80, units = "mm", device = cairo_pdf)
