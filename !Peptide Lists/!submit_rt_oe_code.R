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
  
  detection_folder <- params$detection_folder
})


code$readin_id_all <- quote({
  id_all_precursors_df <-
    here(detection_folder, "id_all_precursors.csv") %>%
    read_id_all_precursors() %>%
    replace_na(list(area_ms1 = 0, area_ms2 = 0)) %>%
    group_by(locator_replicate, locator_peptide) %>%
    arrange(desc(area_ms2), desc(area_ms1), .by_group = TRUE) %>%
    # if multiple charge states, only take the rt from the one with most area
    slice(1) %>%
    ungroup()
})

code$gather_irt_rt <- quote({
  irts <- str_subset(read_lines("!irt_master_list.txt"), "^\\s*#", negate = TRUE)
  
  irts_df <-
    id_all_precursors_df %>%
    select(locator_replicate, rep_name_sky, filename, seq_mod, pre_rt) %>%
    mutate(irt_i = match(seq_mod, irts)) %>%
    filter(!is.na(irt_i), !is.na(pre_rt)) %>%
    arrange(locator_replicate, pre_rt, irt_i)

  irts_df %>%
    count(locator_replicate, rep_name_sky, filename) %>%
    filter(n < 7) %>%
    rowwise() %>%
    group_walk(function(row, grouping) {
      # unfortunately, rowwise walk is performed even with empty dataframe so break out in that case
      if (nrow(row) == 0)
        return()
      if (row$n > 1)
        warning("Replicate ", row$rep_name_sky, " (", row$filename, ") has bad iRT coverage with only ", row$n, " iRT peptides.", call. = FALSE)
      else
        stop("Replicate ", row$rep_name_sky, " (", row$filename, ") has insufficient iRT coverage with only ", row$n, " iRT peptides.", call. = FALSE)
    })
})

code$gather_target_rt <- quote({
  rt_df <-
    id_all_precursors_df %>%
    select(locator_replicate, filename, datetime, seq_mod, pre_rt) %>%
    filter(!(seq_mod %in% irts)) %>%
    mutate(date = as.Date(datetime))
  
  peptides_rt_missing <-
    rt_df %>%
    group_by(seq_mod) %>%
    filter(all(is.na(pre_rt))) %>%
    slice(1) %>%
    pull(seq_mod)
  
  if (!is_empty(peptides_rt_missing))
    warning(paste0("RT missing for: ", paste0(peptides_rt_missing, collapse = ", ")))
  
  if (all(is.na(rt_df$pre_rt)))
    stop("There are no RTs to be submitted.")
})

code$calc_irt_specification <- quote({
  calc_irt <- function(rt, irt_rts, irt_is) {
    stopifnot(!is.unsorted(irt_rts))
    
    if (is.na(rt))
      return(NA_character_)
    
    irts_i <- first(which(irt_rts >= rt))
    if (is.na(irts_i))
      irts_i <- length(irt_rts)
    if (irts_i == 1L)
      irts_i <- 2L
    
    rt_pre <- irt_rts[irts_i - 1L]
    rt_post <- irt_rts[irts_i]
    irt_id_pre <- irt_is[irts_i - 1L]
    irt_id_post <- irt_is[irts_i]
    
    rt_factor <- (rt - rt_pre) / (rt_post - rt_pre)
    
    paste0(irt_id_pre, ":", irt_id_post, "::", round(rt_factor, 3))
  }
  
  rt_df <- 
    rt_df %>%
    filter(!is.na(pre_rt)) %>%
    left_join(
      irts_df %>%
        select(locator_replicate, irt_rts = pre_rt, irt_is = irt_i) %>%
        chop(cols = c(irt_rts, irt_is)),
      by = join_by(locator_replicate)
    ) %>%
    rowwise() %>%
    mutate(irt = calc_irt(pre_rt, irt_rts, irt_is)) %>%
    ungroup() %>%
    select(seq_mod, irt, irt_date = date) %>%
    summarise(.by = seq_mod, irt = paste0(irt, collapse = ";"), irt_date = paste0(irt_date, collapse = ";"))
})

code$read_opt_data <- quote({
  opt_data_df <-
    opt_folder %>%
    list.files(full.names = TRUE, pattern = "\\.csv$") %>%
    tibble(
      path = .,
      name = str_sub(basename(path), end = -5),
      opt_data = map(path, read_csv, col_types = cols(
        seq_mod = col_character(),
        irt = col_character(),
        irt_date = col_character(),
        .default = col_character()
      ))
    ) %>%
    select(-path) %>%
    unnest(opt_data)
})

code$merge_new_opt_data <- quote({
  opt_data_new_df <- rows_update(opt_data_df, rt_df, by = "seq_mod", unmatched = "ignore")
})

code$write_out_opt_data <- quote({
  opt_bkp_folder <- file.path(opt_folder, "bkp")
  if (!dir.exists(opt_bkp_folder))
    dir.create(opt_bkp_folder)
  
  cur_path <- file.path(opt_bkp_folder, "!all_peptides.csv")
  bkp_today_path <- file.path(opt_folder, "bkp", paste0("!all_peptides_", Sys.Date(), ".csv"))
  if (file.exists(cur_path) && !file.exists(bkp_today_path))
    stopifnot(file.copy(file.path(opt_folder, "!all_peptides.csv"), bkp_today_path))
  
  opt_data_new_df %>%
    filter(!is.na(name)) %>%
    group_by(name) %>%
    group_walk(~ write_csv(.x, paste0(opt_folder, "/", .y$name, ".csv"), na = ""))
})
