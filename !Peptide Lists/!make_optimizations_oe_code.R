# This is the outsourced code file for a similarly named .Rmd file
code <- list()

code$setup <- quote({
  library(conflicted)
  library(here)
  library(tidyverse)
  conflicts_prefer(
    dplyr::filter,
    .quiet = TRUE
  )
  read_csv <- partial(readr::read_csv, lazy = FALSE)
  read_tsv <- partial(readr::read_tsv, lazy = FALSE)
  cat0 <- partial(cat, sep = "")
})

code$readin_precursors <- quote({
  precursors_df <-
    "precursors" %>%
    list.files(full.names = TRUE, pattern = "\\.csv$") %>%
    tibble(path = .) %>%
    mutate(
      name = str_sub(basename(path), end = -5),
      precursors = map(path, read_csv, col_types = cols_only(
        seq_mod = col_character(),
        pre_charge = col_integer()
      ))
    ) %>%
    select(-path) %>%
    unnest(precursors)
})


code$readin_optimizations <- quote({
  if (!dir.exists(opt_folder))
    dir.create(opt_folder)
  readin_folder <- file.path(opt_folder, "readin")
  if (!dir.exists(readin_folder))
    dir.create(readin_folder)
  opt_bkp_folder <- file.path(opt_folder, "bkp")
  if (!dir.exists(opt_bkp_folder))
    dir.create(opt_bkp_folder)
  
  not_na_verifyer <- negate(is.na)
  numeric_weak_verifyer <- function(x) is.na(x) | x == "weak" | !is.na(suppressWarnings(as.numeric(x)))
  irt_verifyer <- function(x) is.na(x) | map_lgl(str_split(x, pattern = coll(";")), ~ all(str_detect(.x, "^\\d+:\\d+::-?\\d+\\.?\\d*$")))
  
  opt_columns_df <- tribble(
    ~name,                          ~type,                    ~type_blank,                   ~confirm_func,                                               ~na_to, ~readin_group,
    "seq_mod",                      col_character(),          character(),                   not_na_verifyer,                                             NULL,   NA,
    "pre_charge",                   col_integer(),            integer(),                     not_na_verifyer,                                             NULL,   NA,
    "preferred_precursor",          col_factor(c(1:99, "n")), factor(levels = c(1:99, "n")), NULL,                                                        "n",    1,
    "preferred_precursor_faims",    col_factor(c(1:99, "n")), factor(levels = c(1:99, "n")), NULL,                                                        "n",    2,
    "nce_top_i",                    col_character(),          character(),                   numeric_weak_verifyer,                                       NULL,   3,
    "nce_seq_cov",                  col_character(),          character(),                   numeric_weak_verifyer,                                       NULL,   3,
    "nce_rawfile",                  col_character(),          character(),                   NULL,                                                        NULL,   3,
    "rf_lens",                      col_character(),          character(),                   numeric_weak_verifyer,                                       NULL,   4,
    "rf_lens_rawfile",              col_character(),          character(),                   NULL,                                                        NULL,   4,
    "faims_cv",                     col_character(),          character(),                   numeric_weak_verifyer,                                       NULL,   5,
    "faims_cv_rawfile",             col_character(),          character(),                   NULL,                                                        NULL,   5,
    "irt",                          col_character(),          character(),                   irt_verifyer,                                                NULL,   6,
    "irt_date",                     col_character(),          character(),                   ~ is.na(.x) | !is.na(as.Date(.x)),                           NULL,   6,
  )
  
  col_readers <- do.call(cols_only, deframe(select(opt_columns_df, name, type)))
  col_verifyers <- deframe(select(opt_columns_df, name, confirm_func)) %>% discard(is.null) %>% map(as_mapper)
  col_na_out <- deframe(select(opt_columns_df, name, na_to)) %>% discard(is.null) %>% map(~ function(x) replace_na(x, .x))
  
  id_cols <- c("seq_mod", "pre_charge")
  id_col_symbols <- syms(id_cols)
  
  opt_data_df <-
    opt_folder %>%
    list.files(full.names = TRUE, pattern = "\\.csv$") %>%
    tibble(
      path = .,
      name = str_sub(basename(path), end = -5)
    ) %>%
    rowwise(name) %>%
    reframe(read_csv(path, col_types = col_readers)) %>%
    # ensure all columns are present
    bind_rows(tibble(!!!deframe(select(opt_columns_df, name, type_blank))))
  
  iwalk(col_verifyers, ~ if(!all(.x(opt_data_df[[.y]]))) stop("Column verifyer failed: ", .y))
})


code$consolidate_optimizations <- quote({
  if (only_use_all_peptide_table) {
    opt_data_consolidated_df <- select(filter(opt_data_df, name == "!all_peptides"), -name)
  } else {
    errors <- FALSE
    
    val_or_na <- function(x) {
      vals <- unique(x)
      if (is_scalar_atomic(vals))
        return(vals)
      na.omit(vals)
    }
    
    opt_data_consolidated_df <-
      opt_data_df %>%
      group_by(!!!id_col_symbols) %>%
      group_modify(~ {
        summarised_df <- try(reframe(.x, across(-name, val_or_na)))
        if (is(summarised_df, "tryError") || nrow(summarised_df) > 1) {
          print(.x)
          errors <<- TRUE
        }
        
        summarised_df
      }) %>%
      ungroup()
    
    if (errors)
      stop("Inconsistent data found.")
  }
})


code$list_homeless_precursors <- quote({
  homeless_precursors <-
    opt_data_consolidated_df %>%
    anti_join(precursors_df, by = join_by(!!!id_col_symbols)) %>%
    distinct(!!!id_col_symbols) %>%
    chop(pre_charge) %>%
    rowwise() %>%
    mutate(precursors = paste0(seq_mod, "+[", paste0(pre_charge, collapse = ","), "]")) %>%
    pull(precursors)
  
  if (length(homeless_precursors) > 0)
    warning("The following precursors are homeless: ", paste0(homeless_precursors, "; "))
})


code$integrate_readin_data <- quote({
  readin_groups <-
    opt_columns_df %>%
    select(readin_group, name) %>%
    filter(!is.na(readin_group)) %>%
    group_by(readin_group) %>%
    chop(name) %>%
    deframe()

  readin_df_list <-
    readin_folder %>%
    list.files(full.names = TRUE, pattern = "\\.csv$") %>%
    set_names() %>%
    map(read_csv, col_types = col_readers)
  
  walk(readin_df_list, ~ stopifnot(id_cols %in% colnames(.x)))
  
  if (length(readin_df_list) > 0) {
    cat0("Reading in lists:\n")
    for (path in names(readin_df_list))
      cat0("\t", path, "\n")
  }
  
  for (path in names(readin_df_list)) {
    readin_df <- readin_df_list[[path]]
    if (nrow(readin_df) == 0)
      next()
    
    for (readin_group in readin_groups) {
      readin_avail <- readin_group %in% colnames(readin_df)
      if (!all(readin_avail)) {
        if (any(readin_avail))
          stop("Readin group is incomplete: ", paste0(readin_group, collapse = ", "), "; ", path)
        else
          next()
      }
      
      # define symbols for use in the following steps
      readin_group <- syms(readin_group)
      readin_master <- readin_group[[1]]
      readin_master_y <- sym(paste0(readin_master, ".y"))
      
      readin_overwrite_df <-
        readin_df %>%
        select(!!!id_col_symbols, !!!readin_group) %>%
        filter(readin_force_overwrite | !is.na(!!readin_master)) %>%
        left_join(
          select(opt_data_consolidated_df, !!!id_col_symbols, !!readin_master),
          by = join_by(!!!id_col_symbols),
          suffix = c("", ".y"),
          relationship = "one-to-one"
        ) %>%
        # if the new data is weak, only overwrite old data if it is weak, too
        filter(readin_force_overwrite | !!readin_master != "weak" | (!!readin_master == "weak" & (is.na(!!readin_master_y) | !!readin_master_y == "weak"))) %>%
        select(-!!readin_master_y)
      
      opt_data_consolidated_df <- rows_update(opt_data_consolidated_df, readin_overwrite_df, by = id_cols)
    }
  }
})


code$expand_for_write_out <- quote({
  opt_data_new_df <- full_join(precursors_df, opt_data_consolidated_df, by = join_by(!!!id_col_symbols))
  
  # empty fields of some columns are filled with default data
  for (col in names(col_na_out)) {
    opt_data_new_df[[col]] <- col_na_out[[col]](opt_data_new_df[[col]])
  }
  
  opt_data_new_df <- mutate(opt_data_new_df, preferred_precursor_faims = replace(preferred_precursor, "Y" == preferred_precursor, 1))
  opt_data_new_df <- mutate(opt_data_new_df, preferred_precursor = replace(preferred_precursor, "Y" == preferred_precursor, 1))
})


code$write_out_backup <- quote({
  cur_path <- file.path(opt_folder, "!all_peptides.csv")
  bkp_today_path <- file.path(opt_bkp_folder, paste0("!all_peptides_", Sys.Date(), ".csv"))
  if (file.exists(cur_path) && !file.exists(bkp_today_path))
    stopifnot(file.copy(file.path(opt_folder, "!all_peptides.csv"), bkp_today_path))
})


code$write_out_data <- quote({
  write_csv(distinct(select(opt_data_new_df, -name)), file.path(opt_folder, "!all_peptides.csv"), na = "")
  
  opt_data_new_df %>%
    filter(!is.na(name)) %>%
    group_by(name) %>%
    group_walk(~ {
      cat0(.y$name, "\n")
      write_csv(.x, paste0(opt_folder, "/", .y$name, ".csv"), na = "")
    })
  
  stopifnot(file.remove(names(readin_df_list)))
})
