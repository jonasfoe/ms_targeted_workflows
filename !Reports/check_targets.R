
library(conflicted)
suppressPackageStartupMessages(library(tidyverse))
library(MSTargetedWorkflows)
conflicts_prefer(
  dplyr::filter,
  .quiet = TRUE
)

check_failed <- FALSE

id_all_precursors_df <-
  read_id_all_precursors("id_all_precursors.csv") %>%
  filter(!is.na(filename))

target_peptide_count <-
  id_all_precursors_df %>%
  filter(!(protein %in% c("iRTs", "iRT", "PRTC", "prtc_jptrt_26"))) %>%
  distinct(locator_replicate, seq_mod_sil) %>%
  count(locator_replicate) %>%
  deframe()

ms_files <- set_names(unique(id_all_precursors_df$filename))

ms_files_origin <- find_rawfile(ms_files, base_path = "E:/OE Raw data")

gen_header_rds(ms_files_origin, avoid_regeneration = TRUE)

compounds_df <-
  read_extracted_headers_rds(strip_ms_extension(ms_files)) %>%
  filter(ms_level == 2, scan_mode == "Full") %>%
  distinct(filename_noext = filename, compound) %>%
  mutate(compound = seq_from_rawfile_to_unimod(compound)) %>%
  # add all the metadata for the full join
  left_join(distinct(id_all_precursors_df, locator_replicate, rep_name_sky, filename, filename_noext), by = join_by(filename_noext), relationship = "many-to-many")

check_df <-
  id_all_precursors_df %>%
  mutate(compound = paste0(seq_mod_sil, str_dup("+", pre_charge))) %>%
  distinct(locator_replicate, locator_precursor, rep_name_sky, filename, filename_noext, compound) %>%
  mutate(in_skyline = TRUE) %>%
  full_join(mutate(compounds_df, in_method = TRUE), by = join_by(locator_replicate, rep_name_sky, filename, filename_noext, compound), relationship = "one-to-one") %>%
  replace_na(list(in_skyline = FALSE, in_method = FALSE))

check_df %>%
  group_by(locator_replicate, rep_name_sky, filename) %>%
  group_walk(function(data, grouping) {
    cat(grouping$rep_name_sky, " - ", grouping$filename, "\n")
    cat("\t", "Peptides targeted count: ", target_peptide_count[grouping$locator_replicate], "\n")
    
    only_in_method <- filter(data, in_method, !in_skyline)$compound
    only_in_skyline <- filter(data, !in_method, in_skyline)$compound
    
    if (length(only_in_method) > 0) {
      cat("\t", "Some targets in the method are missing in skyline:", "\n")
      walk(only_in_method, ~ cat("\t\t", .x, "\n"))
    }
    
    if (length(only_in_skyline) > 0) {
      cat("\t", "Some targets are in skyline despite not being targeted in the method:", "\n")
      walk(only_in_skyline, ~ cat("\t\t", .x, "\n"))
    }
    
    if (length(only_in_method) == 0 && length(only_in_skyline) == 0)
      cat("\t", "All good!", "\n")
    else
      check_failed <<- TRUE
  })

if (check_failed)
  stop("Issues were found in the check!", call. = FALSE)
