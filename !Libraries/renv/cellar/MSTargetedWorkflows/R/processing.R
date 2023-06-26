get_process_count <- function() {
  core_count <- parallel::detectCores()
  if (core_count > 1)
    core_count <- core_count - 1L
  core_count
}

as_tmp_file <- function(lines, ext = "") {
  file_path <- tempfile(fileext = ext)
  writeLines(lines, file_path)
  file_path
}

#' @export
process_parallel <- function(command, args) {
  p_count <- length(args)
  p_thread_count <- min(p_count, get_process_count())

  print_output <- function() {
    p_list_curr <- list_processes()
    for (p_out_i in seq_along(p_list_curr)) {
      out <- p_list_curr[[p_out_i]]$read_output_lines()
      err <- p_list_curr[[p_out_i]]$read_error_lines()
      if (!is_empty(out))
        cat(paste0("[", p_out_i, "] ", out, "\n"))
      if (!is_empty(err))
        cat(crayon::red(paste0("[", p_out_i, "] ERROR: ", err, "\n")))
    }
  }

  list_processes <- function() keep(p_list, is, "process")

  count_alive <- function() {
    list_processes() %>%
      keep(~ .x$is_alive()) %>%
      length()
  }

  p_list <- set_names(list_along(args), names(args))

  p_i <- 1
  while (p_i <= p_count) {
    p_list[[p_i]] <- processx::process$new(command, args[[p_i]], echo_cmd = TRUE, stdout = "|", stderr = "|", windows_verbatim_args = TRUE)
    # cat(paste0("[", p_i, "] ", capture.output(print(p_list[[p_i]])), "\n"))

    while (count_alive() == p_thread_count) {
      Sys.sleep(0.25)
      print_output()
    }

    p_i <- p_i + 1
  }

  while (count_alive() > 0) {
    Sys.sleep(0.25)
    print_output()
  }

  print_output()

  errors <-
    p_list %>%
    map_dbl(~ .x$get_exit_status()) %>%
    .[. != 0]

  if (!is_empty(errors)) {
    stop(paste0("Process ", names(errors), " had exit status: ", errors, "\n"))
  }
}

#' @export
gen_header_rds <- function(x, dir = ".", avoid_regeneration = FALSE) {
  x_files <- ensure_extension(x, ".raw")
  names(x_files) <- file.path(dir, paste0(str_sub(basename(x_files), end = -5), ".rds"))

  MsRawAccessVersion <- utils::packageVersion("MsRawAccess")
  apply_versioning <- function(x) {
    attr(x, "MsRawAccessVersion") <- MsRawAccessVersion
    x
  }

  if (avoid_regeneration) {
    to_process <- map_lgl(names(x_files), function(rdsname) {
      if (!file.exists(rdsname))
        return(TRUE)
      version <- attr(readRDS(rdsname), "MsRawAccessVersion")
      if (is.null(version) || version != MsRawAccessVersion)
        return(TRUE)
      FALSE
    })

    x_files <- x_files[to_process]
  }

  count <- length(x_files)
  if (count > 0) {
    cl <- parallel::makeCluster(min(get_process_count(), count))
    x_files %>%
      parallel::parLapplyLB(cl = cl, fun = MsRawAccess::extract_scan_headers_full, chunk.size = 1) %>%
      map(apply_versioning) %>%
      iwalk(saveRDS)
    parallel::stopCluster(cl)
  }

  invisible()
}

#' @export
gen_charges_rds <- function(x, dir = ".", avoid_regeneration = FALSE) {
  x_files <- ensure_extension(x, ".raw")
  names(x_files) <- file.path(dir, paste0(str_sub(basename(x_files), end = -5), "_charges", ".rds"))

  MsRawAccessVersion <- utils::packageVersion("MsRawAccess")
  apply_versioning <- function(x) {
    attr(x, "MsRawAccessVersion") <- MsRawAccessVersion
    x
  }

  if (avoid_regeneration) {
    to_process <- map_lgl(names(x_files), function(rdsname) {
      if (!file.exists(rdsname))
        return(TRUE)
      version <- attr(readRDS(rdsname), "MsRawAccessVersion")
      if (is.null(version) || version != MsRawAccessVersion)
        return(TRUE)
      FALSE
    })

    x_files <- x_files[to_process]
  }

  count <- length(x_files)
  if (count > 0) {
    cl <- parallel::makeCluster(min(get_process_count(), count))
    x_files %>%
      parallel::parLapplyLB(cl = cl, fun = MsRawAccess::extract_ms1_charge_states, chunk.size = 1) %>%
      map(apply_versioning) %>%
      iwalk(saveRDS)
    parallel::stopCluster(cl)
  }

  invisible()
}
