
library(conflicted)
suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
library(MSTargetedWorkflows)
library(ragg)
library(svglite)
conflicts_prefer(
  dplyr::filter,
  .quiet = TRUE
)

theme_set(
  theme_classic(base_size = 7.5, base_family = "sans") +
    theme(
      plot.subtitle = element_text(size = 7),
      axis.text     = element_text(size = 7),
      legend.text   = element_text(size = 7),
      strip.text   = element_text(size = 7.5),
    ) +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      legend.background = element_rect(fill = "transparent"), # legend bg
      legend.box.background = element_rect(fill = "transparent", colour = NA), # legend panel bg
      legend.key = element_rect(fill = "transparent", colour = NA) # legend key bg
    )
)

dot_size <- 3.5

dotp_target <- 0.85

ggsave <- partial(ggsave, bg = "transparent", dpi = 300, units = "mm", limitsize = FALSE)

id_all_precursors_df <-
  read_id_all_precursors("id_all_precursors.csv") %>%
  mutate(has_multi_charge = n_distinct(pre_charge) > 1, .by = locator_peptide) %>%
  mutate(has_multi_hl = n_distinct(is_heavy) > 1, .by = locator_precursor_nohl) %>%
  mutate(compound = paste0(seq_mod_3letter, str_dup("+", ifelse(has_multi_charge, pre_charge, 0))))

id_all_transitions_df <-
  read_id_all_transitions("id_all_transitions.csv") %>%
  filter(scan_type == "Full Ms2")

dotp_selfmade_df <-
  id_all_transitions_df %>%
  summarise(.by = c(locator_replicate, locator_precursor), dotp = calc_nsca(transition_i[transition_is_quantitative], transition_i_lib[transition_is_quantitative])) %>%
  replace_na(list(dotp = 0))

id_all_precursors_df <-
  id_all_precursors_df %>%
  select(-dotp) %>%
  left_join(dotp_selfmade_df, by = join_by(locator_replicate, locator_precursor), relationship = "one-to-one") %>%
  mutate(
    detected = dotp >= dotp_target,
    detected_label = factor(detected, levels = TRUE, labels = "detected")
  )

labels_replicate <- deframe(distinct(id_all_precursors_df, locator_replicate, rep_name_sky))
labels_filename <- deframe(distinct(id_all_precursors_df, locator_replicate, filename_noext))
labels_protein <- deframe(distinct(id_all_precursors_df, locator_protein, protein))
locator_precursor_nohl <- deframe(distinct(id_all_precursors_df, locator_precursor_nohl, compound))

dotp_row_count <- length(levels(id_all_precursors_df$locator_replicate))
dotp_col_count <- 34

dotp_plots_repname <-
  id_all_precursors_df %>%
  filter(!has_multi_hl | !is_heavy) %>%
  group_by(locator_protein, protein) %>%
  mutate(n = (unclass(fct_drop(locator_precursor_nohl)) - 1) %/% dotp_col_count) %>%
  group_by(n, .add = TRUE) %>%
  group_map(function(data, grouping) {
    n_count <- length(levels(fct_drop(data$locator_precursor_nohl)))
    n_start <- grouping$n * dotp_col_count + 1
    
    title <- paste0(grouping$protein, " (", n_start, "-", n_start + n_count - 1, ")")

    data %>%
      ggplot(aes(x = locator_precursor_nohl, y = fct_rev(locator_replicate))) +
      geom_tile(aes(fill = dotp), width = 1, height = 1) +
      geom_text(aes(color = detected_label), label = "\u2022", size = dot_size) +
      scale_x_discrete(labels = locator_precursor_nohl, drop = TRUE) +
      scale_y_discrete(labels = labels_replicate, drop = FALSE) +
      scale_color_manual(values = c(detected = "red"), na.translate = FALSE, drop = FALSE) +
      scale_fill_viridis_c(limits = c(0, 1)) +
      coord_fixed(expand = FALSE, xlim = c(0.5, dotp_col_count + 1 + 0.5)) +
      labs(x = NULL, y = title, fill = "NSA", color = "Indication") +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.line = element_blank()
      )
  })

dotp_plots_filename <-
  id_all_precursors_df %>%
  filter(!has_multi_hl | !is_heavy) %>%
  group_by(locator_protein, protein) %>%
  mutate(n = (unclass(fct_drop(locator_precursor_nohl)) - 1) %/% dotp_col_count) %>%
  group_by(n, .add = TRUE) %>%
  group_map(function(data, grouping) {
    n_count <- length(levels(fct_drop(data$locator_precursor_nohl)))
    n_start <- grouping$n * dotp_col_count + 1
    
    title <- paste0(grouping$protein, " (", n_start, "-", n_start + n_count - 1, ")")
    
    data %>%
      ggplot(aes(x = locator_precursor_nohl, y = fct_reorder(locator_replicate, datetime, .desc = TRUE))) +
      geom_tile(aes(fill = dotp), width = 1, height = 1) +
      geom_text(aes(color = detected_label), label = "\u2022", size = dot_size) +
      scale_x_discrete(labels = locator_precursor_nohl, drop = TRUE) +
      scale_y_discrete(labels = labels_filename, drop = FALSE) +
      scale_color_manual(values = c(detected = "red"), na.translate = FALSE, drop = FALSE) +
      scale_fill_viridis_c(limits = c(0, 1)) +
      coord_fixed(expand = FALSE, xlim = c(0.5, dotp_col_count + 1 + 0.5)) +
      labs(x = NULL, y = title, fill = "NSA", color = "Indication") +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.line = element_blank()
      )
  })

dotp_output_col_count <- 120
dotp_plot_count <- length(dotp_plots_repname)

make_design <- function(n) {
  if (n == 0)
    return("1#")
  return(paste0(c("AB", map_chr(LETTERS[seq_len(n - 1) + 2], ~ paste0("A", .x))), collapse = "\n"))
}
design <- make_design(dotp_plot_count)

dotp_plots_repname_merged <-
  wrap_plots(
    dotp_plots_repname %>%
      map(~ .x + coord_fixed(expand = FALSE, xlim = c(0.5, dotp_output_col_count + 0.5))) %>%
      append(list(guide_area()), after = 0),
    guides = "collect",
    design = design,
    heights = c(rep(100, dotp_plot_count)), widths = c(0.1, 1)
  ) &
  theme(legend.position = "left")

dotp_plots_filename_merged <-
  wrap_plots(
    dotp_plots_filename %>%
      map(~ .x + coord_fixed(expand = FALSE, xlim = c(0.5, dotp_output_col_count + 0.5))) %>%
      append(list(guide_area()), after = 0),
    guides = "collect",
    design = design,
    heights = c(rep(100, dotp_plot_count)), widths = c(0.1, 1)
  ) &
  theme(legend.position = "left")

width <- 500
height <- 100 * dotp_plot_count

saveRDS(dotp_plots_repname, "dotp_plots_repname.rds")
saveRDS(dotp_plots_filename, "dotp_plots_filename.rds")

ggsave("dotp_plots_repname.png", dotp_plots_repname_merged, width = width, height = height, device = agg_png)
ggsave("dotp_plots_repname.svg", dotp_plots_repname_merged, width = width, height = height, device = svglite, fix_text_size = FALSE)
ggsave("dotp_plots_filename.png", dotp_plots_filename_merged, width = width, height = height, device = agg_png)
ggsave("dotp_plots_filename.svg", dotp_plots_filename_merged, width = width, height = height, device = svglite, fix_text_size = FALSE)

# ggpreview(dotp_plots_repname_merged, width = width, height = height, device = agg_png, bg = "transparent", dpi = 300, units = "mm", limitsize = FALSE)
