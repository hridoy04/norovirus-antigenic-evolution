#  ──────────────────────────────────────────────────────────────────────────
#  Plots for noro_nodeall_long.rds: predicted vs observed VP1 mutations.
#
#  Everything meant for tweaking lives in the CONFIG block below. The plot
#  functions read from it, so changing a colour or a domain boundary in one
#  place updates every figure.
#
#    d     <- readRDS("~/noro_nodeall_long.rds")
#    muts  <- noro_mutations(d)        # long: one row per mutation
#    hits  <- noro_hits(d)             # matching summary at 3 stringencies
#    plot_hit_bars(d)
#    plot_protein_map(d)
#  ──────────────────────────────────────────────────────────────────────────

library(tidyverse)


# ── CONFIG ────────────────────────────────────────────────────────────────

# Validated 3-slot categorical palette (CVD dE 9.2, normal-vision dE 24.0).
# Swap these for your own and re-check before publishing.
NORO_PAL <- c(
  observed  = "#2a78d6",   # blue
  predicted = "#eb6834",   # orange
  shared    = "#1baf7a"    # aqua
)

# Ink and chrome. Kept off the series colours on purpose.
NORO_INK <- c(
  primary   = "#0b0b0b",
  secondary = "#52514e",
  muted     = "#8a8984",
  grid      = "#e6e5e1",
  domain    = "#f0efec",
  domain_hi = "#e2e1dc"
)

# GII.4 VP1 domain boundaries, inclusive. P2 ends at 406 to agree with
# NORO_REGIONS in noro_mutation.R; the second P1 segment starts at 407.
VP1_LENGTH <- 540
VP1_DOMAINS <- tribble(
  ~domain, ~start, ~end,
  "S",          1,   225,
  "P1",       226,   278,
  "P2",       279,   406,
  "P1",       407,   540
)

# Blockade epitopes A-E and G. Same definition as NORO_SITES in
# noro_mutation.R and VP1_EPITOPES in noro_deepnode.R -- keep them in step.
VP1_EPITOPES <- list(
  A = c(294:298, 368, 372, 373),
  B = c(333, 389),
  C = c(339:341, 375:378),
  D = c(393:398),
  E = c(407:414),
  G = c(352, 355:357, 359, 364)
)


# ── DATA PREP ─────────────────────────────────────────────────────────────

#  Long form: one row per node x mutation, with the substitution parsed out.
#  `kind` is "observed" (mutation_real, anchored to RefSeq) or "predicted"
#  (mutation_pred, anchored to the node sequence).
noro_mutations <- function(d) {
  rx <- "^([A-Z])([0-9]+)([A-Z*])$"

  d %>%
    select(node, id, date, clinical_state,
           observed = mutation_real, predicted = mutation_pred) %>%
    pivot_longer(c(observed, predicted), names_to = "kind", values_to = "mut") %>%
    unnest_longer(mut) %>%
    filter(!is.na(mut)) %>%
    mutate(
      ancestral = str_match(mut, rx)[, 2],
      position  = as.integer(str_match(mut, rx)[, 3]),
      derived   = str_match(mut, rx)[, 4],
      kind      = factor(kind, levels = c("observed", "predicted"))
    )
}

#  The ancestral residue is anchored differently in the two sets, so exact
#  string matching is only reported for reference -- `pos_aa` and `position`
#  are the defensible keys.
NORO_STRINGENCY <- c(
  exact      = "Exact string",
  pos_aa     = "Position + AA",
  position   = "Position only"
)

noro_keys <- function(muts) {
  muts %>%
    mutate(
      exact    = mut,
      pos_aa   = str_c(position, derived),
      position = as.character(position)
    ) %>%
    pivot_longer(all_of(names(NORO_STRINGENCY)),
                 names_to = "stringency", values_to = "key") %>%
    mutate(stringency = factor(NORO_STRINGENCY[stringency],
                               levels = NORO_STRINGENCY))
}

#  Counts of distinct keys in each set and the overlap, per stringency,
#  pooled over nodes (`scope = "pooled"`) and per node (`scope = "node"`).
noro_hits <- function(d) {
  keyed <- noro_mutations(d) %>% noro_keys()

  pooled <- keyed %>%
    group_by(stringency, kind) %>%
    summarise(keys = list(unique(key)), .groups = "drop") %>%
    pivot_wider(names_from = kind, values_from = keys) %>%
    mutate(
      shared    = map2_int(observed, predicted, ~ length(intersect(.x, .y))),
      observed  = lengths(observed),
      predicted = lengths(predicted)
    ) %>%
    select(stringency, observed, predicted, shared)

  per_node <- keyed %>%
    group_by(stringency, node) %>%
    summarise(
      n_shared = length(intersect(key[kind == "observed"], key[kind == "predicted"])),
      .groups = "drop"
    ) %>%
    group_by(stringency) %>%
    summarise(
      nodes_total  = n(),
      nodes_hit    = sum(n_shared > 0),
      hit_events   = sum(n_shared),
      .groups = "drop"
    )

  list(pooled = pooled, per_node = per_node)
}


# ── THEME ─────────────────────────────────────────────────────────────────

theme_noro <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(colour = NORO_INK[["primary"]]),
      plot.title       = element_text(face = "bold", size = rel(1.15),
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(colour = NORO_INK[["secondary"]],
                                      size = rel(0.9), margin = margin(b = 12)),
      plot.caption     = element_text(colour = NORO_INK[["muted"]],
                                      size = rel(0.75), hjust = 0,
                                      margin = margin(t = 10)),
      axis.title       = element_text(colour = NORO_INK[["secondary"]],
                                      size = rel(0.85)),
      axis.text        = element_text(colour = NORO_INK[["secondary"]],
                                      size = rel(0.85)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = NORO_INK[["grid"]], linewidth = 0.3),
      legend.position  = "top",
      legend.justification = "left",
      legend.title     = element_blank(),
      legend.key.size  = unit(10, "pt"),
      legend.text      = element_text(colour = NORO_INK[["secondary"]],
                                      size = rel(0.85)),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}


# ── PLOT 1: hit counts ────────────────────────────────────────────────────

#  Distinct observed keys, predicted keys and their overlap at each
#  stringency. `metric = "nodes"` instead counts nodes with >=1 hit.
plot_hit_bars <- function(d,
                          metric = c("keys", "nodes"),
                          pal    = NORO_PAL,
                          title  = NULL) {
  metric <- match.arg(metric)
  h <- noro_hits(d)

  if (metric == "keys") {
    dat <- h$pooled %>%
      pivot_longer(c(observed, predicted, shared),
                   names_to = "set", values_to = "n")
    ttl <- title %||% "Observed and predicted VP1 substitutions overlap weakly"
    sub <- "Distinct substitutions pooled across all 129 nodes, by matching stringency"
    ylb <- "Distinct substitutions"
  } else {
    dat <- h$per_node %>%
      transmute(stringency,
                observed  = nodes_total,
                shared    = nodes_hit) %>%
      pivot_longer(c(observed, shared), names_to = "set", values_to = "n")
    ttl <- title %||% "Nodes with at least one matching prediction"
    sub <- "Out of 129 deep nodes, by matching stringency"
    ylb <- "Nodes"
  }

  dat <- dat %>%
    mutate(set = factor(set, levels = names(pal),
                        labels = c("Observed", "Predicted", "Shared")[
                          match(names(pal), names(NORO_PAL))]) %>% fct_drop())

  ggplot(dat, aes(stringency, n, fill = set)) +
    # 2px-equivalent gap between adjacent bars
    geom_col(position = position_dodge(width = 0.78), width = 0.68) +
    geom_text(aes(label = n, group = set),
              position = position_dodge(width = 0.78),
              vjust = -0.6, size = 3.1,
              colour = NORO_INK[["secondary"]]) +
    scale_fill_manual(values = unname(pal)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = ttl, subtitle = sub, x = NULL, y = ylb) +
    theme_noro() +
    theme(panel.grid.major.x = element_blank())
}


# ── PLOT 2: protein map ───────────────────────────────────────────────────

#  VP1 as a 1-540 axis: domain bands behind, observed substitutions as
#  lollipops above the baseline and predicted below, stem height = number of
#  nodes carrying that position. Positions found in both sets are coloured
#  with the `shared` slot and labelled.
plot_protein_map <- function(d,
                             match_on = c("pos_aa", "position", "exact"),
                             label_shared = TRUE,
                             epitopes = VP1_EPITOPES,
                             pal      = NORO_PAL,
                             title    = NULL) {
  match_on <- match.arg(match_on)
  muts <- noro_mutations(d)

  key_of <- function(x) switch(match_on,
    exact    = x$mut,
    pos_aa   = str_c(x$position, x$derived),
    position = as.character(x$position)
  )
  muts <- muts %>% mutate(key = key_of(.))

  shared_keys <- intersect(
    muts$key[muts$kind == "observed"],
    muts$key[muts$kind == "predicted"]
  )

  # one row per kind x position: stem height = nodes carrying it
  stems <- muts %>%
    group_by(kind, position) %>%
    summarise(n_nodes = n_distinct(node),
              is_shared = any(key %in% shared_keys),
              muts = str_c(sort(unique(mut)), collapse = ", "),
              .groups = "drop") %>%
    mutate(y = if_else(kind == "observed", n_nodes, -n_nodes))

  ymax <- max(stems$n_nodes)
  top  <- ymax * 1.30           # headroom for the domain labels

  # P2 shaded; S and P1 left on the surface. Drawn full height so the bands
  # never compete with the marks for the baseline.
  dom <- VP1_DOMAINS %>%
    mutate(fill = if_else(domain == "P2",
                          NORO_INK[["domain_hi"]], NORO_INK[["domain"]]),
           mid  = (start + end) / 2)

  epi <- tibble(epitope = names(epitopes), position = epitopes) %>%
    unnest_longer(position)

  ggplot() +
    # full-height domain bands
    geom_rect(data = dom,
              aes(xmin = start, xmax = end, ymin = -top, ymax = top),
              fill = dom$fill) +
    geom_text(data = dom, aes(mid, top * 0.94, label = domain),
              size = 3, fontface = "bold", colour = NORO_INK[["secondary"]]) +
    geom_hline(yintercept = 0, colour = NORO_INK[["muted"]], linewidth = 0.3) +
    # epitope residues, in a clear lane under the domain labels
    geom_segment(data = epi,
                 aes(x = position, xend = position,
                     y = top * 0.79, yend = top * 0.87),
                 linewidth = 0.45, colour = NORO_INK[["secondary"]]) +
    annotate("text", x = 2, y = top * 0.83, label = "epitopes",
             hjust = 0, size = 2.7, colour = NORO_INK[["muted"]]) +
    # stems
    geom_segment(data = stems,
                 aes(x = position, xend = position, y = 0, yend = y,
                     colour = interaction(kind, is_shared)),
                 linewidth = 0.4) +
    geom_point(data = stems,
               aes(position, y, colour = interaction(kind, is_shared)),
               size = 1.15) +
    scale_colour_manual(
      values = c(
        "observed.FALSE"  = pal[["observed"]],
        "predicted.FALSE" = pal[["predicted"]],
        "observed.TRUE"   = pal[["shared"]],
        "predicted.TRUE"  = pal[["shared"]]
      ),
      breaks = c("observed.FALSE", "predicted.FALSE", "observed.TRUE"),
      labels = c("Observed", "Predicted", "Shared"),
      guide  = guide_legend(override.aes = list(linewidth = 1.2, size = 2))
    ) +
    # shared labels, observed side only to avoid doubling
    {if (label_shared)
      ggrepel::geom_text_repel(
        data = filter(stems, is_shared, kind == "observed"),
        aes(position, y, label = position),
        size = 2.4, colour = NORO_INK[["primary"]],
        segment.colour = NORO_INK[["muted"]], segment.size = 0.2,
        min.segment.length = 0, box.padding = 0.25,
        direction = "y", nudge_y = ymax * 0.1, max.overlaps = Inf
      )
    } +
    annotate("text", x = 2, y = c(top * 0.62, -top * 0.62),
             label = c("observed ↑", "predicted ↓"),
             hjust = 0, size = 3, colour = NORO_INK[["muted"]]) +
    scale_x_continuous(limits = c(0, VP1_LENGTH),
                       breaks = seq(0, VP1_LENGTH, 60),
                       expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(-top, top), labels = abs,
                       expand = expansion(mult = 0)) +
    labs(
      title = title %||% "Predicted substitutions rarely land on observed sites",
      subtitle = str_glue(
        "VP1 codon position; stem height = nodes carrying the substitution. ",
        "Matched on {match_on}."),
      x = "VP1 codon position", y = "Nodes",
      caption = "Shaded band: P2 domain. Ticks: blockade epitope residues."
    ) +
    theme_noro() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = NORO_INK[["grid"]],
                                            linewidth = 0.3))
}
