

library(tidyverse)
library(ggtree)
library(tidytree)
library(TDbook)
library(ggnewscale)
library(ape)
library(treeio)
library(phangorn)


find_deep_nodes <- function(tree) {
  
  library(treeio)
  library(dplyr)
  library(tidyr)
  
  # ── 1. Unnest the MCC tree tibble ──────────────────────────────────────────
  tree_tbl <- tree %>%
    as_tibble() %>%
    mutate(
      location.set      = lapply(location.set, as.character),
      location.set.prob = lapply(location.set.prob, as.numeric)
    ) %>%
    unnest(c(location.set, location.set.prob))
  
  # ── helper: given a set of node IDs, find the deep ancestor per group ──────
  # A "deep ancestor" is the node in a cluster with NO ancestors also in that
  # cluster — i.e. the root of that clade.
  get_deep_ancestors <- function(nodes) {
    if (length(nodes) == 0) return(integer(0))
    
    # Group consecutive node IDs (phylo nodes are numbered in traversal order,
    # so clades tend to be consecutive runs)
    groups <- split(nodes, cumsum(c(1, diff(nodes) != 1)))
    
    deep <- lapply(groups, function(grp) {
      # For every node in the group, find its ancestors that are also in the group
      ancestors_in_grp <- lapply(grp, function(nd) {
        anc <- ancestor(tree, nd)
        anc[anc %in% grp]
      })
      # The deep ancestor is the one with NO in-group ancestors
      grp[lengths(ancestors_in_grp) == 0]
    })
    
    unlist(deep, use.names = FALSE)
  }
  
  # ── 2. Internal nodes only ─────────────────────────────────────────────────
  internal_tbl <- tree_tbl %>% filter(is.na(label))
  
  # ── 3. Pure acute: prob == 1 for "acute" ──────────────────────────────────
  acute_nodes <- internal_tbl %>%
    filter(location.set == "acute", location.set.prob == 1) %>%
    pull(node)
  
  # ── 4. Pure chronic: prob == 1 for "chronic" ──────────────────────────────
  chronic_nodes <- internal_tbl %>%
    filter(location.set == "chronic", location.set.prob == 1) %>%
    pull(node)
  
  # ── 5. Mixed: any node where no single location has prob == 1 ─────────────
  #    i.e. the node appears in the unnested table but never with prob 1
  nodes_with_certain <- internal_tbl %>%
    filter(location.set.prob == 1) %>%
    pull(node) %>%
    unique()
  
  all_internal_nodes <- internal_tbl %>% pull(node) %>% unique()
  
  mixed_nodes <- setdiff(all_internal_nodes, nodes_with_certain)
  
  # ── 6. Find deep ancestors for each category ──────────────────────────────
  list(
    acute   = get_deep_ancestors(acute_nodes),
    chronic = get_deep_ancestors(chronic_nodes),
    mixed   = get_deep_ancestors(mixed_nodes)
  )
}



noro_mcc_beast <-  read.beast('noro7.tree')


result <- find_deep_nodes(noro_mcc_beast)

result$acute    # deep ancestor nodes for pure-acute clades
result$chronic  # deep ancestor nodes for pure-chronic clades
result$mixed    # deep ancestor nodes for mixed/uncertain clades

ggtree(noro_mcc_beast) +
  #geom_tiplab(aes(color = location)) + 
  geom_nodelab(aes(label = location)) +
  geom_hilight(mapping=aes(subset = node %in% result$chronic), fill='darkgreen', type="rect")


#make nice time-calibrated tree with colored clades



ggtree(noro_new_mcc, ladderize = T, mrsd="2024-01-07") + scale_x_continuous(breaks=seq(1970, 2025, 5)) 
noro_new_tree1_subtype <- ggtree(noro_new_mcc, ladderize = T, mrsd="2024-01-07") %<+% noro_new_mcc_metadata + 
  #geom_nodepoint(aes(color = location, subset = !isTip)) + 
  geom_tippoint(aes(color = capsid.subtype),  size = 0.5) + 
  scale_x_continuous(breaks=seq(1970, 2025, 5)) + 
  labs(caption="Years") + 
  scale_color_manual(name = "Subtype/Variant", values = subtype_pal,
                     na.value = "#cccccc", guide = guide_legend(ncol = 1)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_tree2() +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 18),
    legend.title     = element_text(size = 20, face = "bold"),
    # legend.key.size = 4,
    #plot.title       = element_text(size = 12, face = "bold"),
    #axis.text.x      = element_text(size = 7)
  )
noro_new_tree1_subtype + geom_nodelab(aes(label = 2024.0164 - as.numeric(height), subset = !isTip), size = 1, color = "black", hjust = -0.3)
noro_new_tree1_subtype <- noro_new_tree1_subtype + new_scale_fill() 
noro_new_tree1_subtype + geom_nodepoint(aes(color = location, subset = !isTip)) 

n4t_new <- gheatmap(noro_new_tree1_subtype, noro_new_mcc_metadata %>% select(clinical_trait), 
                    offset=.01, width=.05, colnames =F, color =NA) + 
  scale_fill_manual(name = "Clinical Trait", values = trait_pal, 
                    na.value = "black") +
  #scale_x_ggtree() +
  theme_tree2() +
  # theme(
  #   legend.position  = "right",
  #   legend.text      = element_text(size = 18),
  #   legend.title     = element_text(size = 20, face = "bold"),
  #   plot.title       = element_text(size = 24, face = "bold"),
  #   axis.text.x      = element_blank(), 
  # ) +
  ggtitle("Norovirus GII.4 VP1 phylogeny with clinical traits and subtypes (n = 1183)")

ggsave('norotree_new_v2.pdf', n4t_new, width = 18.38, height = 13.94,  units  = "in",
       device = cairo_pdf)


#get nodes and offspring details

# Tabulate every internal node with its date, clinical state and descendant tips.
# `mrsd` is the most recent sampling date, as "YYYY-MM-DD" or a decimal year.

get_node_offspring <- function(tree, mrsd) {
  
  mrsd_dec <- if (is.numeric(mrsd)) mrsd else lubridate::decimal_date(as.Date(mrsd))
  
  tree_tbl  <- as_tibble(tree)
  tip_names <- tree_tbl %>% select(node, id = label)
  
  tree_tbl %>%
    filter(is.na(label)) %>%
    transmute(
      node,
      clinical_state = location,
      date           = mrsd_dec - as.numeric(height),
      offspring      = map(node, ~ offspring(tree, .x, type = "tips"))  # tips only, internal nodes excluded
    ) %>%
    unnest_longer(offspring) %>%
    left_join(tip_names, by = c("offspring" = "node"))
}

noro_new_mcc <- read.beast('/Users/rubayetalam/noro7.tree')

nndf <- get_node_offspring(noro_new_mcc, mrsd = "2024-01-07")
nndf %>% write_csv('nndf.csv')

nndf %>% filter(node %in% result$chronic)
length(unique(nndf$id))


#get antigenicity pred

cna <- read_csv('cna.csv')

cna_all <- cna %>% group_by(virus_name_protein_name, position...1) %>% 
  slice_max(significant_antigenicity_median, na_rm = T) %>% 
  filter(label...4 == label...10) %>% 
  select(virus_name_protein_name, position...1, ref...6, alt...7, significant_antigenicity_median, label...4, SurfaceAccessibility, aa_change_likelihood) %>% 
  rename(
    virus_name = virus_name_protein_name,
    position = position...1,
    epal = significant_antigenicity_median,
    label = label...4,
    ref = ref...6,
    alt = alt...7
    ) %>% 
  filter(SurfaceAccessibility >= 0.15)

write_csv(cna_all, 'cna_all.csv')

cna_all

noro_deepnode_impmut <- cna_all %>%
  group_by(virus_name) %>%
  summarise(pos = list(position), label = list(label))
ndnim <- noro_deepnode_impmut %>% rename(mut = label)

write_csv(noro_deepnode_impmut, 'noro_deepnode_impmut.csv')


#make a single df

#mutation_df

nm1183 <- read_csv('noro_mutation_1183.csv')
nm1183


noro_deepnode_all_df <- nm1183 %>% right_join(nndfa, by = 'id') %>% 
  group_by(node, ancestor_node) %>% 
  mutate(id_list = list(id)) %>% ungroup() %>% 
  distinct(ancestor_node, .keep_all = T) %>% 
  left_join(ndnim, by = c('ancestor_node' = 'virus_name')) %>% 
  filter(lengths(pos) > 0) %>% 
  rename(mutation_real = mutation, 
         mutation_pred = mut) 

write(noro_deepnode_all_df, 'noro_deepnode_all_df.csv')

noro_deepnode_all_df_rds <- noro_deepnode_all_df %>% mutate(mutation_real = str_extract_all(mutation_real, "[A-Z][0-9]+[A-Z]"))

saveRDS(noro_deepnode_all_df_rds, "noro_nodeall_long.rds")

noro_deepnode_all_df %>% select(ancestor_node, clinical_state) %>% write_csv('deepnode.csv')





