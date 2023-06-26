
library(conflicted)
suppressPackageStartupMessages(library(tidyverse))
library(MSTargetedWorkflows)
library(ggpp)
conflicts_prefer(
  dplyr::filter,
  ggpp::annotate,
  .quiet = TRUE,
)

fontsize <- 5
width <- 150
height <- 20
units <- "mm"

id_all_precursors_df <-
  read_id_all_precursors("id_all_precursors.csv")

id_all_transtions_df <-
  id_all_precursors_df %>%
  left_join(read_id_all_transitions("id_all_transitions.csv"), by = join_by(locator_replicate, locator_precursor), relationship = "one-to-many")

pdf(NULL)
coverage_plots_df <-
  id_all_transtions_df %>%
  filter(!is.na(pre_rt)) %>%
  group_by(locator_replicate, locator_precursor, rep_name_sky, protein, seq_mod_3letter, pre_charge) %>%
  filter(transition_is_coeluting, transition_i > 0, .preserve = TRUE) %>%
  group_modify(function(data, grouping) {
    cov_plot <- plot_peptide_sequence_coverage(
      grouping$seq_mod_3letter, grouping$pre_charge,
      data$transition_type, data$transition_ordinal,
      width = width, height = height, units = "mm", vjust = 1, fontsize = fontsize,
      mod_fontsize_factor = 0.8
    )

    cov_plot_repname_annotated <- cov_plot + annotate("text_npc", label = paste0(grouping$rep_name_sky, " - ", grouping$protein), npcx = 0.01, npcy = 0.95, hjust = 0, vjust = 1, size = fontsize * 0.6)

    tibble(cov_plot = list(cov_plot), cov_plot_repname_annotated = list(cov_plot_repname_annotated))
  }) %>%
  ungroup()
invisible(dev.off())

list(
  width = 210,
  height = 297,
  units = "mm",
  cov_plots_df = coverage_plots_df
) %>%
  saveRDS("sequence_coverage_plots.rds")

pdf("sequence_coverage_plots.pdf", width = width / 25.4, height = height / 25.4, bg = "transparent")
walk(coverage_plots_df$cov_plot_repname_annotated, print)
invisible(dev.off())

