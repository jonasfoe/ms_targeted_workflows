# This is the outsourced code file for a similarly named .Rmd file
code <- list()

code$setup <- quote({
  library(conflicted)
  library(here)
  library(tidyverse)
  library(MSTargetedWorkflows)
  library(patchwork)
  conflicts_prefer(
    dplyr::filter,
    .quiet = TRUE
  )
  read_csv <- partial(readr::read_csv, lazy = FALSE)
  read_tsv <- partial(readr::read_tsv, lazy = FALSE)
  cat0 <- partial(cat, sep = "")
  top_i <- sym("nce_top_i")
  seq_cov <- sym("nce_seq_cov")
  suppressWarnings(rm("seq_include_only", inherits = FALSE))
  suppressWarnings(rm("seq_include_extra", inherits = FALSE))
})

code$process_settings <- quote({
  stopifnot(is.symbol(nce_type), as.character(nce_type) %in% c("nce_top_i", "nce_seq_cov"))
  nce_label <- str_sub(as.character(nce_type), start = 5)
})

code$readin_irts <- quote({
  irts_ref <- irt_list$prtc_jptrt_26
  
  irt_fallback_df <-
    here(irt_fallback_folder, "id_all_precursors.csv") %>%
    read_id_all_precursors() %>%
    select(seq_mod, rt_fallback = pre_rt) %>%
    filter(seq_mod %in% irts_ref)
  
  stopifnot(check_df_single_cases(irt_fallback_df, single_case_cols = everything(), by = seq_mod))
  stopifnot(all(irts_ref %in% irt_fallback_df$seq_mod))
  
  irt_latest_df <-
    here(irt_latest_folder, "id_all_precursors.csv") %>%
    read_id_all_precursors() %>%
    select(seq_mod, rt_latest = pre_rt) %>%
    filter(seq_mod %in% irts_ref)
  
  stopifnot(check_df_single_cases(irt_latest_df, single_case_cols = everything(), by = seq_mod))
  stopifnot(nrow(irt_latest_df) >= 2)
  
  irt_df <-
    irt_fallback_df %>%
    left_join(irt_latest_df, by = join_by(seq_mod), relationship = "one-to-one") %>%
    mutate(rt_fitted = rt_latest) %>%
    mutate(seq_mod = factor(seq_mod, levels = irts_ref)) %>%
    arrange(seq_mod)
})

code$irt_fallback_merge <- quote({
  rt_available <- which(!is.na(irt_df$rt_latest))
  
  for (i in 1:nrow(irt_df)) {
    if (i %in% rt_available)
      next()
    i_rt_pre <- last(rt_available[rt_available < i])
    i_rt_post <- first(rt_available[rt_available > i])
    
    if (is.na(i_rt_pre)) {
      i_rt_pre <- i_rt_post
      i_rt_post <- first(rt_available[rt_available > i_rt_pre])
    } else if (is.na(i_rt_post)) {
      i_rt_post <- i_rt_pre
      i_rt_pre <- first(rt_available[rt_available < i_rt_post])
    }
    
    stopifnot(!is.na(i_rt_pre), !is.na(i_rt_post))
    
    rt_fallback_pre <- irt_df$rt_fallback[i_rt_pre]
    rt_fallback_post <- irt_df$rt_fallback[i_rt_post]
    rt_fallback <- irt_df$rt_fallback[i]
    
    rt_delta <- (rt_fallback - rt_fallback_pre) / (rt_fallback_post - rt_fallback_pre)
    irt_df$rt_fitted[i] <- irt_df$rt_latest[i_rt_pre] + (irt_df$rt_latest[i_rt_post] - irt_df$rt_latest[i_rt_pre]) * rt_delta
  }
})

code$compute_peptide_set_name <- quote({
  peptide_set_name <- peptide_set
  if (length(peptide_set_name) == 0)
    peptide_set_name <- str_trunc(paste0(as.character(seq_include_extra), collapse = "_"), width = 42, side = "center", ellipsis = "__")
  if (length(peptide_set_name) > 1)
    peptide_set_name <- str_trunc(paste0(peptide_set_name, collapse = "_"), width = 42, side = "center", ellipsis = "__")
})

code$readin_precursors <- quote({
  precursors_irts_df <-
    here("!Peptide Lists", "precursors", "prtc_jptrt_26.csv") %>%
    read_csv(col_types = cols_only(
      n = col_double(),
      list_number = col_double(),
      seq = col_character(),
      seq_sil = col_character(),
      seq_mod = col_character(),
      seq_mod_sil = col_character(),
      seq_sil_masses = col_character(),
      seq_mod_masses = col_character(),
      seq_mod_sil_masses = col_character(),
      formula_heavy = col_character(),
      pre_charge = col_double(),
      mz_heavy = col_double()
    )) %>%
    mutate(
      across(c(seq, seq_sil, seq_mod, seq_mod_sil, seq_sil_masses, seq_mod_masses, seq_mod_sil_masses), fct_inorder),
      n = unclass(seq_mod)
    )
  
  read_precursors <- function(path) {
    read_csv(path, col_types = cols_only(
      n = col_integer(),
      list_number = col_integer(),
      seq = col_character(),
      seq_sil = col_character(),
      seq_mod = col_character(),
      seq_mod_sil = col_character(),
      seq_sil_masses = col_character(),
      seq_mod_masses = col_character(),
      seq_mod_sil_masses = col_character(),
      formula_heavy = col_character(),
      formula_light = col_character(),
      pre_charge = col_integer(),
      mz_heavy = col_double(),
      mz_light = col_double()
    ))
  }
  
  precursors <-
    peptide_set %>%
    map_chr(~ here("!Peptide Lists", "precursors", paste0(.x, ".csv"))) %>%
    map(read_precursors)
  
  precursors_extra_df <- imap(
    seq_include_extra,
    function(seq, set)
      filter(read_precursors(here("!Peptide Lists", "precursors", paste0(set, ".csv"))), seq == !!seq)
  )
  
  precursors_df <-
    bind_rows(!!!precursors, !!!precursors_extra_df) %>%
    # only the first heavy labelling is kept
    distinct(seq_mod, pre_charge, .keep_all = TRUE) %>%
    mutate(
      across(c(seq, seq_sil, seq_mod, seq_mod_sil, seq_sil_masses, seq_mod_masses, seq_mod_sil_masses), fct_inorder),
      n = unclass(seq_mod)
    )
})

code$readin_precursor_optimizations <- quote({
  weak_to_na <- function(x) replace(x, x == "weak", NA_character_)
  
  optimizations_irts_df <-
    here("!Peptide Lists", "optimized_oe", "prtc_jptrt_26.csv") %>%
    read_csv(col_types = cols_only(
      seq_mod = col_character(),
      pre_charge = col_double(),
      preferred_precursor = col_factor(c(1:99, "n")),
      nce_top_i = col_character(),
      nce_seq_cov = col_character(),
      rf_lens = col_character(),
      faims_cv = col_character()
    )) %>%
    mutate(across(preferred_precursor, ~ as.numeric(replace(as.character(.x), .x == "n", NA_character_)))) %>%
    mutate(across(c(nce_top_i, nce_seq_cov, rf_lens), weak_to_na)) %>%
    mutate(across(c(nce_top_i, nce_seq_cov, rf_lens), as.numeric)) %>%
    mutate(nce = !!nce_type, .before = nce_top_i)
  
  optimizations_df <-
    here("!Peptide Lists", "optimized_oe", "!all_peptides.csv") %>%
    read_csv(col_types = cols_only(
      seq_mod = col_character(),
      pre_charge = col_double(),
      preferred_precursor = col_factor(c(1:99, "n")),
      preferred_precursor_faims = col_factor(c(1:99, "n")),
      nce_top_i = col_character(),
      nce_seq_cov = col_character(),
      rf_lens = col_character(),
      faims_cv = col_character(),
      irt = col_character(),
      irt_date = col_character()
    )) %>%
    mutate(across(c(preferred_precursor, preferred_precursor_faims), ~ as.numeric(replace(as.character(.x), .x == "n", NA_character_)))) %>%
    mutate(across(c(nce_top_i, nce_seq_cov, rf_lens), weak_to_na)) %>%
    mutate(across(c(nce_top_i, nce_seq_cov, rf_lens), as.numeric)) %>%
    mutate(nce = !!nce_type, .before = nce_top_i)
})

code$create_schedules <- quote({
  irt_to_rt <- function(irt) {
    stopifnot(!is.na(irt))
    
    irt_rts <- irt_df$rt_fitted
    
    irt_split <- str_match(irt, "^(\\d+):(\\d+)::(-?\\d+\\.?\\d*)$")
    irt_i_pre <- as.integer(irt_split[, 2])
    irt_i_post <- as.integer(irt_split[, 3])
    irt_delta <- as.numeric(irt_split[, 4])
    
    rt_pre <- irt_rts[irt_i_pre]
    rt_post <- irt_rts[irt_i_post]
    
    rt_pre + (rt_post - rt_pre) * irt_delta
  }
  
  apply_schedule_filters <- function(df) {
    for (fun in schedule_filters) df <- as_mapper(fun)(df)
    df
  }
  
  schedule_irts_df <-
    precursors_irts_df %>%
    left_join(optimizations_irts_df, by = join_by(seq_mod, pre_charge), relationship = "one-to-one") %>%
    filter(seq_mod %in% irts, preferred_precursor == 1) %>%
    replace_na(list(nce = nce_default, rf_lens = rf_default, faims_cv = faims_default)) %>%
    left_join(select(irt_df, seq_mod, rt = rt_fitted), by = join_by(seq_mod), relationship = "one-to-one")
  
  stopifnot(all(irts %in% schedule_irts_df$seq))
  
  schedule_df <-
    precursors_df %>%
    left_join(optimizations_df, by = join_by(seq_mod, pre_charge), relationship = "many-to-one") %>%
    filter(
      !is.na(irt),
      !preferred_precursors_only |
        (!faims_use & !is.na(preferred_precursor) & preferred_precursor <= preferred_precursors_include_priority) |
        (faims_use & !is.na(preferred_precursor_faims) & preferred_precursor_faims <= preferred_precursors_include_priority),
      !nce_use | (!nce_only_optimized | !is.na(nce)),
      !rf_use  | (!rf_only_optimized | !is.na(rf_lens)),
      !faims_use | (!faims_only_optimized | !is.na(faims_cv))
    ) %>%
    replace_na(list(nce = nce_default, rf_lens = rf_default, faims_cv = faims_default)) %>%
    mutate(irt = str_split(irt, ";")) %>%
    unchop(irt) %>%
    mutate(rt = map_dbl(irt, irt_to_rt))
  
  precursors_lost <- discard(unique(precursors_df$seq_mod), `%in%`, schedule_df$seq_mod)
  
  schedule_df <- apply_schedule_filters(schedule_df)
  if (exists("seq_include_only", inherits = FALSE) && is.character(seq_include_only) && length(seq_include_only) > 0)
    schedule_df <- filter(schedule_df, seq %in% seq_include_only)
  
  if (!is_empty(precursors_lost))
    warning("Peptides lost due to missing optimization: ", paste0(precursors_lost, collapse = ", "))
})

code$quad_isolation_rules <- quote({
  oe_quadrupole_iso_df <- tribble(
    ~mz, ~iso_w_mz,
    400,       0.4,
    700,       0.7,
    1000,         1,
    1500,       1.5,
    2000,         2,
    2500,         3,
  )
  
  # get optimal precursor isolation window according to thermo OE manual
  get_isolation_window <- function(mz) {
    stopifnot(mz <= 2500)
    
    i <- which.max(oe_quadrupole_iso_df$mz >= mz)
    
    oe_quadrupole_iso_df$iso_w_mz[i]
  }
})

code$create_dynrt_schedule <- quote({
  schedule_irts_dynrt_out_df <-
    schedule_irts_df %>%
    filter(seq_mod %in% irt_list$pierce_15) %>%
    mutate(rt_window = if_else(seq == "SSAAPPPPPR", rt_window_prtc_first, rt_window_prtc)) %>%
    select(
      `Peptide Name` = seq_mod,
      `m/z` = mz_heavy,
      `RT Time (min)` = rt,
      `RT Window (min)` = rt_window
    )
  
  write_csv(schedule_irts_dynrt_out_df, file.path(out_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_dynRT.csv")))
})

code$create_irt_schedule <- quote({
  schedule_irts_out_df <-
    schedule_irts_df %>%
    mutate(
      Compound = paste0(seq_mod_sil, str_dup("+", pre_charge)),
      Adduct = "+H",
      iso_w_mz = map_dbl(mz_heavy, get_isolation_window),
      rt_window = if_else(seq == "SSAAPPPPPR", rt_window_prtc_first, rt_window_irt)
    ) %>%
    select(
      Compound,
      Formula = formula_heavy,
      Adduct,
      `m/z` = mz_heavy,
      z = pre_charge,
      `RT Time (min)` = rt,
      `Window (min)` = rt_window,
      `Isolation Window (m/z)` = iso_w_mz,
      `HCD Collision Energy (%)` = nce,
      `FAIMS CV (V)` = faims_cv
    )
  
  if (!faims_use)
    schedule_irts_out_df <- select(schedule_irts_out_df, -`FAIMS CV (V)`)
  
  write_csv(schedule_irts_out_df, file.path(out_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_irt_", nce_label, "_schedule.csv")))
})

code$create_schedule <- quote({
  irt_coverage_min <- min(schedule_irts_df$rt) - rt_window / 2
  irt_coverage_max <- max(schedule_irts_df$rt) + rt_window / 2
  rt_window_outside_diff <- rt_window_outside_irt - rt_window
  rt_window_transition_border_lower <- irt_coverage_min - rt_window_outside_diff
  rt_window_transition_border_upper <- irt_coverage_min + rt_window_outside_diff
  
  calc_rt_window <- function(rt) {
    if (between(rt, irt_coverage_min, irt_coverage_max))
      return(rt_window)
    
    if (rt < rt_window_transition_border_lower || rt > rt_window_transition_border_upper)
      return(rt_window_outside_irt)
    
    if (rt < irt_coverage_min)
      return(rt_window + (irt_coverage_min - rt))
    
    # if (rt > irt_coverage_max)
    return(rt_window + rt - irt_coverage_max)
  }
    
  merge_rt_overlaps <- function(df) {
    df %>%
      group_by(across(c(-rt_start, -rt_end))) %>%
      group_modify(function(rt_df, info) {
        if (nrow(rt_df) == 1)
          return(rt_df)
        rt_df <- arrange(rt_df, rt_end)
        rt_gap <- head(rt_df$rt_end, -1) - tail(rt_df$rt_start, -1)
        has_overlap <- which(rt_gap >= 0)
        if (is_empty(has_overlap))
          return(rt_df)
        tibble(
          rt_start = rt_df$rt_start[-(has_overlap + 1)],
          rt_end = rt_df$rt_end[-has_overlap]
        )
      }) %>%
      ungroup()
  }
  
  schedule_out_proto_df <-
    schedule_df %>%
    mutate(
      rt_window = map_dbl(rt, calc_rt_window),
      rt_start = pmax(rt - rt_window / 2, 0),
      rt_end = rt + (rt_window / 2)
    ) %>%
    select(-irt, -rt, -rt_window) %>%
    merge_rt_overlaps()
  
  schedule_out_light_df <- 
    schedule_out_proto_df %>%
    mutate(
      Compound = paste0(seq_mod, str_dup("+", pre_charge)),
      Adduct = "+H",
      iso_w_mz = map_dbl(mz_light, get_isolation_window),
    ) %>%
    select(
      Compound,
      Formula = formula_light,
      Adduct,
      `m/z` = mz_light,
      z = pre_charge,
      `t start (min)` = rt_start,
      `t stop (min)` = rt_end,
      `Isolation Window (m/z)` = iso_w_mz,
      `HCD Collision Energy (%)` = nce,
      `RF Lens (%)` = rf_lens,
      `FAIMS CV (V)` = faims_cv
    )
  
  if (!nce_use)
    schedule_out_light_df <- select(schedule_out_light_df, -`HCD Collision Energy (%)`)
  if (!rf_use)
    schedule_out_light_df <- select(schedule_out_light_df, -`RF Lens (%)`)
  if (!faims_use)
    schedule_out_light_df <- select(schedule_out_light_df, -`FAIMS CV (V)`)
  
  write_csv(schedule_out_light_df, file.path(out_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_light_", nce_label, "_schedule.csv")))
  
  schedule_out_heavy_df <-
    schedule_out_proto_df %>%
    mutate(
      Compound = paste0(seq_mod_sil, str_dup("+", pre_charge)),
      Adduct = "+H",
      iso_w_mz = map_dbl(mz_heavy, get_isolation_window),
    ) %>%
    select(
      Compound,
      Formula = formula_heavy,
      Adduct,
      `m/z` = mz_heavy,
      z = pre_charge,
      `t start (min)` = rt_start,
      `t stop (min)` = rt_end,
      `Isolation Window (m/z)` = iso_w_mz,
      `HCD Collision Energy (%)` = nce,
      `RF Lens (%)` = rf_lens,
      `FAIMS CV (V)` = faims_cv
    )
  
  if (!nce_use)
    schedule_out_heavy_df <- select(schedule_out_heavy_df, -`HCD Collision Energy (%)`)
  if (!rf_use)
    schedule_out_heavy_df <- select(schedule_out_heavy_df, -`RF Lens (%)`)
  if (!faims_use)
    schedule_out_heavy_df <- select(schedule_out_heavy_df, -`FAIMS CV (V)`)
  
  write_csv(schedule_out_heavy_df, file.path(out_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_heavy_", nce_label, "_schedule.csv")))
  
  schedule_out_light_df %>%
    pull(Compound) %>%
    as.character() %>%
    unique() %>%
    writeLines(file.path(out_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_compounds.txt")))
  
  schedule_out_proto_df %>%
    pull(seq_mod_sil_masses) %>%
    as.character() %>%
    unique() %>%
    writeLines(file.path(out_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_seq_for_skyline.txt")))
})

code$create_schedule_overlap_plot <- quote({
  overlaps_targets_prior_df <- tibble(rt = numeric(), val = numeric())
  if (exists("overlaps_targets_df"))
    overlaps_targets_prior_df <- bind_rows(overlaps_targets_prior_df, overlaps_targets_df)
  
  overlaps_targets_df <-
    bind_rows(
      tibble(rt = 0, val = 0),
      tibble(rt = schedule_out_light_df$`t start (min)`, val = 1),
      tibble(rt = schedule_out_light_df$`t stop (min)`, val = -1),
    )
  
  overlaps_df <- bind_rows(
    irts = bind_rows(
      tibble(rt = 0, val = 0),
      tibble(rt = schedule_irts_out_df$`RT Time (min)` - schedule_irts_out_df$`Window (min)` / 2, val = 1),
      tibble(rt = schedule_irts_out_df$`RT Time (min)` + schedule_irts_out_df$`Window (min)` / 2, val = -1),
    ),
    targets_prior = overlaps_targets_prior_df,
    targets = overlaps_targets_df,
    .id = "set"
  ) %>%
    mutate(set = fct_inorder(set)) %>%
    arrange(rt) %>%
    mutate(.by = set, count = cumsum(val))
  
  overlaps_plot <-
    overlaps_df %>%
    ggplot(aes(rt, count, color = set, linetype = set)) +
    geom_step(direction = "hv") +
    scale_x_continuous(breaks = scales::breaks_extended(Q = 1:5, n = 15)) +
    scale_y_continuous(breaks = 1:1e4, minor_breaks = NULL) +
    scale_color_manual(values = c(irts = "blue", targets_prior = "darkgrey", targets = "black"), breaks = rev(levels(overlaps_df$set))) +
    scale_linetype_manual(values = c(irts = 2, targets_prior = 1, targets = 1), breaks = rev(levels(overlaps_df$set))) +
    labs(title = "MS2 RT Overlaps", subtitle = "Use Ctrl+Alt+R to recalculate based on current settings.", x = "Retention Time (min)", y = "Simultaneous Scheduled Targets (#)")
  
  rt_windows_plot <-
    schedule_out_light_df %>%
    mutate(rt_mean = `t start (min)` + (`t stop (min)` - `t start (min)`) / 2) %>%
    arrange(desc(rt_mean), desc(z)) %>%
    mutate(Compound = fct_inorder(Compound)) %>%
    ggplot(aes(x = `t start (min)`, xend = `t stop (min)`, y = Compound, yend = Compound)) +
    geom_segment(linewidth = 1) +
    scale_x_continuous(breaks = scales::breaks_extended(Q = 1:5, n = 15)) +
    coord_cartesian(xlim = range(overlaps_df$rt)) +
    labs(x = "Retention Time (min)", y = "Compounds")
  
  # workaround for a bug with patchwork
  schedue_overlaps_plot <- function() patchwork::wrap_plots(overlaps_plot, rt_windows_plot, ncol = 1, heights = c(3, 2))
})


code$create_rt_age_plot <- quote({
  plot_rt_age <-
    schedule_df %>%
    mutate(irt_date = str_split(irt_date, coll(";"))) %>%
    unchop(irt_date) %>%
    distinct(seq_mod, irt_date) %>%
    ggplot(aes(y = seq_mod, x = difftime(Sys.Date(), irt_date, unit = "days"))) +
    geom_col(position = position_dodge2(padding = 0.03, preserve = "total")) +
    scale_x_continuous(breaks = scales::breaks_extended(Q = 1:5, n = 15), minor_breaks = NULL) +
    labs(title = "Peptide RT reference age", y = NULL, x = "Age (days)")
})

# code$create_tsim_schedule <- quote({
#   calc_tsim_rt_start <- function(rt) map_dbl(rt - (rt_window_tsim / 2), max, 0)
#   calc_tsim_rt_end <- function(rt) rt + (rt_window_tsim / 2)
#   
#   schedule_tsim_out_light_df <-
#     schedule_df %>%
#     mutate(rt_start = calc_tsim_rt_start(rt), rt_end = calc_tsim_rt_end(rt)) %>%
#     select(-irt, -rt) %>%
#     merge_rt_overlaps() %>%
#     mutate(
#       Compound = paste0(seq_semi_mod_sky, str_dup("+", pre_charge)),
#       Adduct = "+H",
#       # isolation window is based on including M+2
#       iso_w_mz = map_dbl(mz_light + 2 / pre_charge, get_isolation_window) + 4 / pre_charge
#     ) %>%
#     select(
#       Compound,
#       Formula = formula_light,
#       Adduct,
#       `m/z` = mz_light,
#       z = pre_charge,
#       `t start (min)` = rt_start,
#       `t stop (min)` = rt_end,
#       `Isolation Window (m/z)` = iso_w_mz,
#       `RF Lens (%)` = rf_lens,
#       `FAIMS CV (V)` = faims_cv
#     )
#   
#   if (!rf_use)
#     schedule_tsim_out_light_df <- select(schedule_tsim_out_light_df, -`RF Lens (%)`)
#   if (!faims_use)
#     schedule_tsim_out_light_df <- select(schedule_tsim_out_light_df, -`FAIMS CV (V)`)
#   
#   write_csv(schedule_tsim_out_light_df, paste0("schedules/", format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_tsim_light_schedule.csv"))
#   
#   schedule_tsim_out_heavy_df <-
#     schedule_df %>%
#     mutate(rt_start = calc_tsim_rt_start(rt), rt_end = calc_tsim_rt_end(rt)) %>%
#     select(-irt, -rt) %>%
#     merge_rt_overlaps() %>%
#     mutate(
#       Compound = paste0(seq_mod_sky, str_dup("+", pre_charge)),
#       Adduct = "+H",
#       # isolation window is based on including M+2
#       iso_w_mz = map_dbl(mz_heavy + 2 / pre_charge, get_isolation_window) + 4 / pre_charge
#     ) %>%
#     select(
#       Compound,
#       Formula = formula_heavy,
#       Adduct,
#       `m/z` = mz_heavy,
#       z = pre_charge,
#       `t start (min)` = rt_start,
#       `t stop (min)` = rt_end,
#       `Isolation Window (m/z)` = iso_w_mz,
#       `RF Lens (%)` = rf_lens,
#       `FAIMS CV (V)` = faims_cv
#     )
#   
#   if (!rf_use)
#     schedule_tsim_out_heavy_df <- select(schedule_tsim_out_heavy_df, -`RF Lens (%)`)
#   if (!faims_use)
#     schedule_tsim_out_heavy_df <- select(schedule_tsim_out_heavy_df, -`FAIMS CV (V)`)
#   
#   write_csv(schedule_tsim_out_heavy_df, paste0("schedules/", format(Sys.Date(), "%Y%m%d"), "_", peptide_set_name, "_tsim_heavy_schedule.csv"))
# })

# code$create_tsim_schedule_overlap_plot <- quote({
#   plot_tsim_overlaps <-
#     bind_rows(
#       tibble(rt = 0, val = 0),
#       tibble(rt = schedule_tsim_out_light_df$`t start (min)`, val = 1),
#       tibble(rt = schedule_tsim_out_light_df$`t stop (min)`, val = -1)
#     ) %>%
#     arrange(rt) %>%
#     mutate(count = cumsum(val)) %>%
#     ggplot(aes(rt, count)) +
#     geom_step(direction = "hv") +
#     scale_x_continuous(n.breaks = 15) +
#     scale_y_continuous(breaks = 1:9999, minor_breaks = NULL) +
#     ggtitle("tSIM RT Overlaps")
# })

code$calculate_peptide_mixes <- quote({
  mixes_df <-
    readLines(here("!Peptide Lists", "!peptide_mixes.txt")) %>%
    str_subset("^\\s*#", negate = TRUE) %>%
    set_names() %>%
    map_chr(~ here("!Peptide Lists", "precursors", paste0(.x, ".csv"))) %>%
    map(read_csv, col_types = cols_only(
      seq_mod_sil = col_character()
    )) %>%
    bind_rows(.id = "mix") %>%
    distinct()
  
  # hacked together method to get a min amount of mixes for all peptides
  # bad performance but could be ok for now as most peptides are in only one mix
  mixes_mashup_df <-
    schedule_df %>%
    distinct(seq_mod_sil) %>%
    left_join(mixes_df, by = join_by(seq_mod_sil), relationship = "one-to-many")
  
  if (anyNA(mixes_mashup_df$mix)) {
    pep_not_in_mixes_df <- mixes_mashup_df %>% filter(is.na(mix)) %>% distinct(seq_mod_sil)
    warning("Some peptides are not in any mixes (as defined in !Peptide Lists/!peptide_mixes.txt): ", paste0(pep_not_in_mixes_df$seq_mod_sil, collapse = ", "))
    mixes_mashup_df <- anti_join(mixes_mashup_df, pep_not_in_mixes_df, by = join_by(seq_mod_sil))
    cat0("\n")
  }
  
  mixes_forced_df <-
    mixes_mashup_df %>%
    chop(mix) %>%
    filter(lengths(mix) == 1) %>%
    unchop(mix) %>%
    distinct(mix)
  
  mixes_additional <-
    mixes_mashup_df %>%
    anti_join(mixes_forced_df, by = join_by(mix)) %>%
    # remove all peptides that are already covered by the forced mixes
    anti_join(distinct(inner_join(mixes_mashup_df, mixes_forced_df, by = join_by(mix)), seq_mod_sil), by = join_by(seq_mod_sil)) %>%
    chop(mix) %>%
    pull(mix) %>%
    unique() %>%
    expand_grid(!!!.) %>%
    pmap(\(...) sort(unique(c(...)))) %>%
    unique()
  
  if (length(mixes_additional) == 0)
    cat0("1: ", paste0(mixes_forced_df$mix, collapse = ", "), "\n")
  else
    mixes_additional %>%
      map(append, mixes_forced_df$mix) %>%
      map(sort) %>%
      map_chr(paste0, collapse = ", ") %>%
      iwalk(~ cat0(.y, ": ", .x, "\n"))
})
