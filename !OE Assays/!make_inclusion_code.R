# This is the outsourced code file for a similarly named .Rmd file
code <- list()

code$setup <- quote({
  library(conflicted)
  library(here)
  library(tidyverse)
  library(MSTargetedWorkflows)
  conflicts_prefer(
    dplyr::filter,
    .quiet = TRUE
  )
  read_csv <- partial(readr::read_csv, lazy = FALSE)
  read_tsv <- partial(readr::read_tsv, lazy = FALSE)
  cat0 <- partial(cat, sep = "")
  top_i <- sym("nce_top_i")
  seq_cov <- sym("nce_seq_cov")
})

code$process_settings <- quote({
  stopifnot(is.symbol(nce_type), as.character(nce_type) %in% c("nce_top_i", "nce_seq_cov"))
})

code$compute_peptide_set_name <- quote({
  peptide_set_name <- peptide_set
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
  
  precursors_df <-
    peptide_set %>%
    map_chr(~ here("!Peptide Lists", "precursors", paste0(.x, ".csv"))) %>%
    map(read_csv, col_types = cols_only(
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
    )) %>%
    bind_rows() %>%
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
      irt_date = col_character()
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
      nce_top_i = col_character(),
      nce_seq_cov = col_character(),
      rf_lens = col_character(),
      irt = col_character(),
      irt_date = col_character()
    )) %>%
    mutate(across(c(preferred_precursor), ~ as.numeric(replace(as.character(.x), .x == "n", NA_character_)))) %>%
    mutate(across(c(nce_top_i, nce_seq_cov, rf_lens), weak_to_na)) %>%
    mutate(across(c(nce_top_i, nce_seq_cov, rf_lens), as.numeric)) %>%
    mutate(nce = !!nce_type, .before = nce_top_i)
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

code$make_inclusion <- quote({
  apply_inclusion_filters <- function(df) {
    for (fun in inclusion_filters) df <- as_mapper(fun)(df)
    df
  }
  
  inclusion_irts_df <-
    precursors_irts_df %>%
    left_join(optimizations_irts_df, by = join_by(seq_mod, pre_charge), relationship = "many-to-one") %>%
    filter(seq_mod %in% irts, preferred_precursor == 1) %>%
    replace_na(list(nce = nce_default, rf_lens = rf_default))
  
  inclusion_df <-
    precursors_df %>%
    left_join(optimizations_df, by = join_by(seq_mod, pre_charge), relationship = "many-to-one") %>%
    filter(
      !preferred_precursors_only | preferred_precursor <= preferred_precursors_include_priority,
      !remove_discovered_peptides | is.na(irt),
      !nce_only_optimized | !is.na(nce_top_i),
      !rf_only_optimized | !is.na(rf_lens),
    ) %>%
    replace_na(list(nce = nce_default, rf_lens = rf_default))
  
  precursors_lost <- discard(unique(precursors_df$seq_mod), `%in%`, inclusion_df$seq_mod)
  
  inclusion_df <- apply_inclusion_filters(inclusion_df)
  
  if (!is_empty(precursors_lost))
    warning("Peptides lost due to missing optimization: ", paste0(precursors_lost, collapse = ", "))
  
  inclusion_full_df <-
    bind_rows(
      irts = inclusion_irts_df,
      targets = inclusion_df,
      .id = "set"
    ) %>%
    mutate(set = fct_inorder(set)) %>%
    arrange(set, n) %>%
    select(set, seq_mod_sil, seq_mod_sil_masses, pre_charge, formula_heavy, mz_heavy, nce) %>%
    mutate(
      Compound = fct_inorder(paste0(seq_mod_sil, str_dup("+", pre_charge))),
      mz_iso_window = map_dbl(mz_heavy, get_isolation_window)
    )
})

code$write_out <- quote({
  inclusion_out_df <-
    inclusion_full_df %>%
    select(
      Compound,
      Formula = formula_heavy,
      `m/z` = mz_heavy,
      z = pre_charge,
      # `Isolation Window (m/z)` = mz_iso_window,
      `HCD Collision Energy (%)` = nce
    ) %>%
    mutate(Adduct = "+H", .after = Formula)
  
  cat0("Writing out info under name '", peptide_set_name, "'", "\n")
  inclusion_out_name <- paste0(peptide_set_name, "_iDDA_inclusion_", ms1_range_label, ".csv")
  write_csv(inclusion_out_df, file.path(out_dir, inclusion_out_name))
  cat0("Inclusion list: ", inclusion_out_name, "\n")
  
  # inclusion_tsim_out_df <-
  #   inclusion_full_df %>%
  #   mutate(mz_iso_window = map2_dbl(mz_iso_window, 5 / pre_charge, max)) %>%
  #   select(
  #     Compound,
  #     `m/z` = mz_heavy,
  #     z = pre_charge,
  #     `Isolation Window (m/z)` = mz_iso_window
  #   ) %>%
  #   mutate(Formula = "", Adduct = "", .after = Compound)
  # 
  # write_csv(inclusion_tsim_out_df, file.path(out_dir, paste0(peptide_set_name, "_tSIM_inclusion.csv")))
  
  inclusion_full_df %>%
    filter(set == "targets") %>%
    pull(Compound) %>%
    as.character() %>%
    unique() %>%
    writeLines(file.path(out_dir, paste0(peptide_set_name, "_compounds.txt")))
  
  inclusion_full_df %>%
    filter(set == "targets") %>%
    pull(seq_mod_sil_masses) %>%
    as.character() %>%
    unique() %>%
    writeLines(file.path(out_dir, paste0(peptide_set_name, "_seq_for_skyline.txt")))
})

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
    inclusion_df %>%
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
