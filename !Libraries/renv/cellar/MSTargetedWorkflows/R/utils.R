#' @export
str_prepend <- function(x, y) {
  if (is_empty(x)) return(x)
  set_names(str_c(y, x), names(x))
}

#' @export
str_append <- function(x, y) {
  if (is_empty(x)) return(x)
  set_names(str_c(x, y), names(x))
}

#' @export
single_row_df_to_list <- function(df) {
  assert_that(nrow(df) == 1)
  map(df, first)
}

#' @export
swap_names_and_content <- function(x) set_names(names(x), unname(x))

#' @export
ensure_names <- function(...) {
  x <- c(...)
  nm <- names(x) %||% x
  set_names(x, replace(nm, nm == "", x[nm == ""]))
}

#' @export
list_along <- function(x) {
  vector("list", length = length(x))
}

# Check that cols in the `single_case_cols` argument contain only one case each, per grouping as specified in `by`
#' @export
check_df_single_cases <- function(df, single_case_cols = everything(), by = c()) {
  distinct_df <- distinct(df, pick({{ by }}), pick({{ single_case_cols }}))
  counts_df <- count(distinct_df, pick({{ by }}))
  all(pull(counts_df, last_col()) == 1)
}

# Check that all rows are distinct, including only those listed in the `cols` argument
#' @export
check_df_all_distinct <- function(df, cols = everything()) {
  counts_df <- count(df, pick({{ cols }}))
  all(pull(counts_df, last_col()) == 1)
}

# integrate irregular timeseries using trapezoidal rule
# allows for exact limits (off timesteps)
#' @export
integrate_timeseries <- function(x, y, x_lims = c(NA_real_, NA_real_)) {
  assert_that(!is.unsorted(x), length(x) >= 1L, length(y) >= 1L, length(x) == length(y))
  assert_that(is.numeric(x_lims), length(x_lims) == 2L)

  do_start <- FALSE
  do_end <- FALSE

  i_first <- 1L
  i_last <- length(x)

  if (!is.na(x_lims[1])) {
    assert_that(x_lims[1] >= first(x))
    if (x_lims[1] > first(x)) {
      do_start <- TRUE
      i_first <- match(TRUE, x > x_lims[1])
    }
  }
  if (!is.na(x_lims[2])) {
    assert_that(x_lims[2] <= last(x))
    if (!is.na(x_lims[1]))
      assert_that(x_lims[1] <= x_lims[2])
    if (x_lims[2] < last(x)) {
      do_end <- TRUE
      i_last <- match(TRUE, x > x_lims[2]) - 1L
    }
  }

  x_body <- x[i_first:i_last]
  y_body <- y[i_first:i_last]

  x_diff <- diff(x_body)
  x_weights <- c(x_diff, 0) + c(0, x_diff)

  x_result <- sum(x_weights * y_body) / 2

  if (do_start) {
    x1 <- x[i_first - 1L]; x2 <- x[i_first]; y1 <- y[i_first - 1L]; y2 <- y[i_first]
    dydx <- (y2 - y1) / (x2 - x1)
    dx <- x2 - x_lims[1]
    x_result <- x_result + y2 * dx - dydx * dx ^ 2 / 2
  }
  if (do_end) {
    x1 <- x[i_last]; x2 <- x[i_last + 1L]; y1 <- y[i_last]; y2 <- y[i_last + 1L]
    dydx <- (y2 - y1) / (x2 - x1)
    dx <- x_lims[2] - x1
    x_result <- x_result + y1 * dx + dydx * dx ^ 2 / 2
  }
  x_result
}

# normalized spectral contrast angle as per Toprak et al
#' @export
calc_nsca <- function(x, y) {
  x <- sqrt(x) / sqrt(sum(x))
  y <- sqrt(y) / sqrt(sum(y))
  1 - (2 * acos(drop(x %*% y))) / pi
}

#' @export
ensure_extension <- function(x, ext) {
  if (str_sub(ext, 1, 1) != ".")
    ext <- str_c(".", ext)

  has_corr_ext <- str_sub(x, -str_length(ext), -1) == ext

  replace(x, !has_corr_ext, str_append(x[!has_corr_ext], ext))
}

#' @export
strip_ms_extension <- function(x) {
  map_chr(x, ~ {
    if (is.na(.x))
      NA_character_
    else if (str_ends(.x, coll(".raw")))
      str_sub(.x, end = -5)
    else if (str_ends(.x, coll(".mzml")) || str_ends(.x, coll(".mzML")))
      str_sub(.x, end = -6)
    else
      stop("no ms extension to strip")
  })
}

#' @export
find_rawfile <- function(x, base_path) {
  names_bkp <- names(x)

  find_file <- function(filename, date) {
    if (file.exists(filename))
      return(filename)

    filepath <- file.path(base_path, filename)
    if (file.exists(filepath))
      return(filepath)

    filepath <- file.path(base_path, date[1], paste0(date[1], date[2]), paste0(date[2], date[3]), filename)
    if (file.exists(filepath))
      return(filepath)

    return(NA_character_)
  }

  out <-
    tibble(base = str_remove(x, "\\.raw$")) %>%
    mutate(
      filename = str_append(base, ".raw"),
      date = str_match(base, "(\\d\\d\\d\\d)(\\d\\d)(\\d\\d)_(\\d\\d\\d)(?>_[\\w\\d-]+)?$")[, 2:5, drop = FALSE]
    ) %>%
    rowwise() %>%
    summarise(path = find_file(filename, date), .groups = "drop") %>%
    pull(path)

  names(out) <- names_bkp

  out
}

#' @export
irt_list <-
  list(
    pierce_15 = c(
      "SSAAPPPPPR",
      "GISNEGQNASIK",
      "HVLTSIGEK",
      "DIPVPKPK",
      "IGDYAGIK",
      "TASEFDSAIAQDK",
      "SAAGAFGPELSR",
      "ELGQSGVDTYLQTK",
      "GLILVGGYGTR",
      "GILFVGSGVSGGEEGAR",
      "SFANQPLEVVYSK",
      "LTILEELR",
      "NGFILDGFPR",
      "ELASGLSFPVGFK",
      "LSSEAPALFQFDLK"
    ),
    biognosis_11 = c(
      "LGGNEQVTR",
      "GAGSSEPVTGLDAK",
      "VEATFGVDESNAK",
      "YILAGVENSK",
      "TPVISGGPYEYR",
      "TPVITGAPYEYR",
      "DGLDAASYYAPVR",
      "ADVTPADFSEWSK",
      "GTFIIDPGGVIR",
      "GTFIIDPAAVIR",
      "LFLQFGAQGSPFLK"
    ),
    prtc_jptrt_26 = c(
      "YSAHEEHHYDK",
      "HEHISSDYAGK",
      "SSAAPPPPPR",
      "GISNEGQNASIK",
      "HVLTSIGEK",
      "TASGVGGFSTK",
      "DIPVPKPK",
      "IGDYAGIK",
      "TLIAYDDSSTK",
      "TASEFDSAIAQDK",
      "SAAGAFGPELSR",
      "ELGQSGVDTYLQTK",
      "HLTGLTFDTYK",
      "HDTVFGSYLYK",
      "GLILVGGYGTR",
      "SFANQPLEVVYSK",
      "GILFVGSGVSGGEEGAR",
      "ASDLLSGYYIK",
      "TSIDSFIDSYK",
      "LTILEELR",
      "NGFILDGFPR",
      "ELASGLSFPVGFK",
      "LSSEAPALFQFDLK",
      "GDFTFFIDTFK",
      "SLFFIIDGFVK",
      "SLIFFLSTLLK"
    )
  )

#' @export
is_irt <- function(x) {
  x %>%
    map_lgl(~ {
      matches <- map_dbl(irt_list, match, x = .x)
      sum(!is.na(matches)) > 0
    })
}

#' @export
which_irt <- function(x) {
  x %>%
    map_dbl(~ {
      matches <- map_dbl(irt_list, match, x = .x)
      which(!is.na(matches))[1] %||% NA
    }) %>%
    {names(irt_list)[.]}
}

#' @export
get_irt_accession <- function(x) {
  x %>%
    which_irt() %>%
    str_c("irt_", ., "_") %>%
    replace_na("")
}

#' @export
irt_validator <- c(
  pierce_15 = "^[A-Z]{7,16}[RK]\\(Label:13C\\(6\\)15N\\([42]\\)\\)$",
  biognosis_11 = "^[A-Z]{9,14}$"
)

#' @export
has_irt_accession <- function(x) {
  str_detect(x, '^irt_\\w+_\\d+_')
}

#' @export
which_irt_accession <- function(x) {
  str_match(x, '^irt_(\\w+_\\d+)_')[, 2]
}

#' @export
is_invalid_irt <- function(name, seq) {
  is_irt <- has_irt_accession(name)
  ret <- is_irt

  irt_name <- which_irt_accession(name[is_irt])

  ret[is_irt] <- !str_detect(seq[is_irt], irt_validator[irt_name])

  ret
}

skyline_full_fixup <- c(
  "M\\[Oxidation \\(M\\)[^\\]]*\\]" = "M[Oxidation (M)]",
  "C\\[Carbamidomethyl \\(C\\)[^\\]]*\\]" = "C[Carbamidomethyl (C)]"
)

rawfile_compound_to_unimod <- c(
  "M[+15.994915]" = "M[Oxidation]",
  "C[+57.021464]" = "C[Carbamidomethyl]",
  "K[+8.014199]" = "K[Label:13C(6)15N(2)]",
  "R[+10.008269]" = "R[Label:13C(6)15N(4)]",
  "I[+7.017164]" = "I[Label:13C(6)15N(1)]",
  "L[+7.017164]" = "L[Label:13C(6)15N(1)]",
  "P[+6.013809]" = "P[Label:13C(5)15N(1)]",
  "V[+6.013809]" = "V[Label:13C(5)15N(1)]",
  "F[+10.027228]" = "F[Label:13C(9)15N(1)]",
  "A[+4.007099]" = "A[Label:13C(3)15N(1)]"
)

unimod_to_3letter <- c(
  "M[Oxidation]" = "M[Oxi]",
  "C[Carbamidomethyl]" = "C[CAM]",
  "K[Label:13C(6)15N(2)]" = "K[+08]",
  "R[Label:13C(6)15N(4)]" = "R[+10]",
  "I[Label:13C(6)15N(1)]" = "I[+07]",
  "L[Label:13C(6)15N(1)]" = "L[+07]",
  "P[Label:13C(5)15N(1)]" = "P[+06]",
  "V[Label:13C(5)15N(1)]" = "V[+06]",
  "F[Label:13C(9)15N(1)]" = "F[+10]",
  "A[Label:13C(3)15N(1)]" = "A[+04]"
)

#' @export
seq_skyline_full_fixup <- function(seq) {
  str_replace_all(seq, pattern = regex(skyline_full_fixup))
}

#' @export
seq_from_skyline_full_to_unimod <- function(seq) {
  seq %>%
    seq_skyline_full_fixup() %>%
    str_replace_all("([A-Z])\\[([^\\]]+)\\s\\(\\1\\)]", "\\1[\\2]")
}

#' @export
seq_from_rawfile_to_unimod <- function(seq) {
  str_replace_all(seq, pattern = coll(rawfile_compound_to_unimod))
}

#' @export
seq_from_unimod_to_3letter <- function(seq) {
  str_replace_all(seq, pattern = coll(unimod_to_3letter))
}

#' @export
seq_from_skyline_full_to_skyline_3letter <- function(seq) {
  seq %>%
    seq_from_skyline_full_to_unimod() %>%
    str_replace_all(pattern = coll(unimod_to_3letter))
}

#' @export
seq_from_3letter_to_unimod <- function(seq) {
  str_replace_all(seq, pattern = coll(swap_names_and_content(unimod_to_3letter)))
}

#' @export
seq_demod_unimod <- function(x) {
  str_replace_all(x, "\\[[^\\]]*\\]", "")
}

#' @export
sample_df <- function(df, col = everything(), facets = c(), factor = 20) {
  sample_df <-
    df %>%
    group_by(pick({{ facets }})) %>%
    distinct(pick({{ col }})) %>%
    arrange(pick({{ col }})) %>%
    filter(row_number() %% factor == 1)

  semi_join(df, sample_df, by = names(sample_df))
}

#' @export
read_extracted_headers_rds <- function(x, minimal = FALSE) {
  ms_lvl_to_string <- function(ms_lvl) replace(paste0("Ms", ms_lvl),  ms_lvl == 1, "Ms")

  columns = tribble(
    ~pre,                      ~post,                    ~convert_to,
    "origin",                  "origin",                 NA,
    "filename",                "filename",               NA,
    "ScanNumber",              "scan_number",            NA,
    "ScanParent",              "scan_parent_number",     NA,
    "StartTime",               "rt",                     NA,
    "Ion Injection Time (ms)", "it",                     "numeric",
    "HCD Energy",              "nce",                    "numeric",
    "TIC",                     "tic",                    NA,
    "MsOrder",                 "ms_level",               NA,
    "ScanMode",                "scan_mode",              NA,
    "Compound",                "compound",               NA,
    "IsDependent",             "is_dd",                  NA,
    "IsMultiplex",             "is_mpx",                 NA,
    "HasLock",                 "has_lock",               NA,
    "PrecursorMasses",         "pre_mz",                 NA,
    "MS2 Isolation Width",     "pre_window",             "numeric",
    "MS2 Isolation Offset",    "pre_window_offset",      "numeric",
    "Max. Ion Time (ms)",      "max_it",                 "numeric",
    "FT Resolution",           "resolution",             "numeric",
    "AGC Target",              "agc_target",             "numeric",
    "Dynamic RT Shift (min)",  "rt_shift",               "numeric",
    "AGC Fill",                "agc_fill",               "numeric",
    "LM m/z-Correction (ppm)", "lm_mz_corr",             "numeric",
    "Number of LM Found",      "lm_found_count",         "numeric",
  )

  headers_df <-
    x %>%
    ensure_names() %>%
    ensure_extension(".rds") %>%
    map(~
      .x %>%
      readRDS() %>%
      mutate(filename = str_sub(.x, end = -5))
    ) %>%
    bind_rows(.id = "origin")

  columns_available <-
    columns$pre %>%
    set_names(columns$post) %>%
    .[. %in% names(headers_df)]

  # reorder and rename the columns
  if (minimal) {
    headers_df <- select(headers_df, !!!columns_available)
  } else {
    headers_df <- select(headers_df, !!!columns_available, everything())
  }

  headers_df <-
    headers_df %>%
    mutate(
      origin = factor(origin, levels = unique(origin)),
      filename = factor(filename, levels = unique(filename))
    )

  column_coversions <-
    columns %>%
    filter(post %in% names(headers_df), !is.na(convert_to)) %>%
    select(post, convert_to) %>%
    tibble::deframe()

  for (k in names(column_coversions))
    headers_df[[k]] <- as.vector(headers_df[[k]], mode = column_coversions[k])

  if (has_name(headers_df, "pre_mz"))
    headers_df <- mutate(headers_df, pre_mz = str_split(pre_mz, ";") %>% map(as.numeric))


  if (has_name(headers_df, "scan_mode"))
    headers_df <- mutate(headers_df, scan_mode = factor(scan_mode, levels = c("Full", "Sim")))

  if (has_name(headers_df, "ms_level")) {
    headers_df <- mutate(headers_df,
                         ms_level = str_match(ms_level, "^Ms(\\d*)$")[ ,2] %>% replace(. == "", "1") %>% as.numeric()
    )
    if (has_name(headers_df, "scan_mode"))
      headers_df <- mutate(headers_df, scan_type = paste(scan_mode, ms_lvl_to_string(ms_level)))
  }

  headers_df
}

#' @export
read_extracted_charges_rds <- function(x) {
  columns = tribble(
    ~pre,            ~post,
    "origin",        "origin",
    "filename",      "filename",
    "scan_number",   "scan_number",
    "rt",            "rt",
    "charge",        "charge",
    "intensity",     "intensity",
  )

  charges_df <-
    x %>%
    ensure_names() %>%
    ensure_extension(".rds") %>%
    map(~
      .x %>%
      readRDS() %>%
      mutate(filename = str_sub(.x, end = -13))) %>%
    bind_rows(.id = "origin")

  columns_names <- set_names(columns$pre, columns$post)

  # reorder and rename the columns
  charges_df <- select(charges_df, !!!columns_names)

  charges_df <-
    charges_df %>%
    mutate(
      origin = factor(origin, levels = unique(origin)),
      filename = factor(filename, levels = unique(filename))
    )

  charges_df
}

#' @export
read_id_all_precursors <- function(path) {
  columns_df = tribble(
    ~pre, ~post, ~type,
    "File Name", "filename", readr::col_character(),
    "Replicate Name", "rep_name_sky", readr::col_character(),
    "Acquired Time", "datetime", readr::col_character(),
    "Replicate Locator", "locator_replicate", readr::col_factor(),
    "Protein Name", "protein", readr::col_character(),
    "Protein Locator", "locator_protein", readr::col_factor(),
    "Peptide Sequence", "seq", readr::col_character(),
    "Peptide Modified Sequence Full Names", "seq_sky_mod", readr::col_character(),
    "Peptide Modified Sequence Three Letter Codes", "seq_mod_3letter", readr::col_character(),
    "Peptide Locator", "locator_peptide", readr::col_factor(),
    "Precursor Charge", "pre_charge", readr::col_integer(),
    "Modified Sequence Full Names", "seq_sky_mod_sil", readr::col_character(),
    "Modified Sequence Three Letter Codes", "seq_mod_sil_3letter", readr::col_character(),
    "Isotope Label Type", "is_heavy", readr::col_character(),
    "Precursor Mz", "pre_mz", readr::col_double(),
    "Precursor Ion Formula", "pre_formula", readr::col_character(),
    "Best Retention Time", "pre_rt", readr::col_double(),
    "Min Start Time", "pre_rt_min", readr::col_double(),
    "Max End Time", "pre_rt_max", readr::col_double(),
    "Predicted Retention Time", "rt_predicted", readr::col_double(),
    "Total Area Ms1", "area_ms1", readr::col_double(),
    "Total Area Fragment", "area_ms2", readr::col_double(),
    "Isotope Dot Product", "idotp", readr::col_double(),
    "Library Dot Product", "dotp", readr::col_double(),
    "Ratio Dot Product", "rdotp", readr::col_double(),
    "Detection Q Value", "q_value", readr::col_double(),
    "Precursor Locator", "locator_precursor", readr::col_factor(),
  )

  columns_df <- mutate(columns_df, pre = str_remove_all(pre, "\\s"))

  df <- readr::read_csv(path, col_types = do.call(readr::cols, set_names(as.list(columns_df$type), columns_df$pre)), na = c("", "NA", "#N/A"), lazy = FALSE)

  unknown_pos <- nrow(columns_df) - 1
  renamer1 <- columns_df %>% slice(1:unknown_pos) %>% filter(pre %in% names(df)) %>% select(post, pre) %>% tibble::deframe()
  renamer2 <- columns_df %>% slice((unknown_pos + 1):n()) %>% filter(pre %in% names(df)) %>% select(post, pre) %>% tibble::deframe()
  df <- select(df, !!!renamer1, everything(), !!!renamer2)

  if (has_name(df, "filename"))
    df <- mutate(df, .after = filename, filename_noext = strip_ms_extension(filename))

  if (has_name(df, "datetime"))
    df <- mutate(df, datetime = as.POSIXct(datetime, format = c("%m/%d/%Y %H:%M:%S")))

  if (has_name(df, "is_heavy"))
    df <- mutate(df, is_heavy = is_heavy == "heavy")

  if (has_name(df, "locator_precursor"))
    df <- mutate(df, .after = pre_charge,
      locator_precursor_nohl = fct_reorder(str_remove(locator_precursor, "\\/(?>light|heavy)(?=\\++$)"), unclass(locator_precursor))
    )

  if (has_name(df, "seq_sky_mod"))
    df <- mutate(df, seq_mod = seq_from_skyline_full_to_unimod(seq_sky_mod), .before = seq_sky_mod)

  if (has_name(df, "seq_sky_mod_sil"))
    df <- mutate(df, seq_mod_sil = seq_from_skyline_full_to_unimod(seq_sky_mod_sil), .before = seq_sky_mod_sil)

  # skyline 3 letter export is often messy so I try to fix some known issues here
  if (has_name(df, "seq_mod_3letter"))
    df <- mutate(df, seq_mod_3letter = seq_from_skyline_full_to_skyline_3letter(seq_mod_3letter))

  if (has_name(df, "seq_mod_sil_3letter"))
    df <- mutate(df, seq_mod_sil_3letter = seq_from_skyline_full_to_skyline_3letter(seq_mod_sil_3letter))

  df
}

#' @export
read_id_all_transitions <- function(path) {
  columns_df = tribble(
    ~pre,                                            ~post,                        ~type,
    "Replicate Locator",                             "locator_replicate",          readr::col_factor(),
    "Precursor Locator",                             "locator_precursor",          readr::col_factor(),
    "Fragment Ion",                                  "transition",                 readr::col_character(),
    "Fragment Ion Type",                             "transition_type",            readr::col_character(),
    "Fragment Ion Ordinal",                          "transition_ordinal",         readr::col_integer(),
    "Product Charge",                                "transition_charge",          readr::col_integer(),
    "Product Mz",                                    "transition_mz",              readr::col_double(),
    "Peak Rank By Level",                            "transition_rank",            readr::col_integer(),
    "Area",                                          "transition_i",               readr::col_double(),
    "Background",                                    "transition_i_bg",            readr::col_double(),
    "Library Rank",                                  "transition_rank_lib",        readr::col_integer(),
    "Library Intensity",                             "transition_i_lib",           readr::col_double(),
    "Quantitative",                                  "transition_is_quantitative", readr::col_logical(),
    "Coeluting",                                     "transition_is_coeluting",    readr::col_logical(),
    "Chromatogram Source",                           "scan_type",                  readr::col_character(),
    "Transition Locator",                            "locator_transition",         readr::col_factor(),
  )

  columns_df <- mutate(columns_df, pre = str_remove_all(pre, "\\s"))

  df <- readr::read_csv(path, col_types = do.call(readr::cols, set_names(as.list(columns_df$type), columns_df$pre)), na = c("", "NA", "#N/A"), lazy = FALSE)

  unknown_pos <- nrow(columns_df) - 1
  renamer1 <- columns_df %>% slice(1:unknown_pos) %>% filter(pre %in% names(df)) %>% select(post, pre) %>% tibble::deframe()
  renamer2 <- columns_df %>% slice((unknown_pos + 1):n()) %>% filter(pre %in% names(df)) %>% select(post, pre) %>% tibble::deframe()
  df <- select(df, !!!renamer1, everything(), !!!renamer2)

  if (has_name(df, "scan_type"))
    df <- mutate(df, scan_type = factor(
      scan_type,
      levels = c("ms1", "sim", "fragment"),
      labels = c("Full Ms", "SIM Ms", "Full Ms2")
    ))

  df <- mutate(df, transition_name = paste0(transition, str_dup("+", transition_charge)), .before = transition)

  df
}

#' @export
read_id_all_chromatograms <- function(path) {
  columns_df = tribble(
    ~pre,                       ~post,                        ~type,
    "Replicate Locator",        "locator_replicate",          readr::col_factor(),
    "Transition Locator",       "locator_transition",         readr::col_factor(),
    "Raw Spectrum Ids",         "scan_number_raw",            readr::col_character(),
    "Raw Times",                "rt_raw",                     readr::col_character(),
    "Raw Intensities",          "i_raw",                      readr::col_character(),
  )

  columns_df <- mutate(columns_df, pre = str_remove_all(pre, "\\s"))

  df <- readr::read_csv(path, col_types = do.call(readr::cols, set_names(as.list(columns_df$type), columns_df$pre)), na = c("", "NA", "#N/A"), lazy = FALSE)

  renamer <- columns_df %>% filter(pre %in% names(df)) %>% select(post, pre) %>% tibble::deframe()
  df <- select(df, !!!renamer, everything())

  df <-
    df %>%
    replace_na(list(scan_number_raw = "", rt_raw = "", i_raw = "")) %>%
    rowwise(c(-scan_number_raw, -rt_raw, -i_raw)) %>%
    reframe(
      scan_number =
        scan_number_raw %>%
        str_split_1(",") %>%
        str_extract("(?<=^0\\.1\\.)\\d+$") %>%
        as.integer,
      rt =
        rt_raw %>%
        str_split_1(",") %>%
        as.numeric,
      i =
        i_raw %>%
        str_split_1(",") %>%
        as.numeric
    ) %>%
    nest(chromatogram_df = c(scan_number, rt, i))

  df
}

#' @export
ggpreview <- function(...) {
  loadNamespace("ggplot2")
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible()
}
