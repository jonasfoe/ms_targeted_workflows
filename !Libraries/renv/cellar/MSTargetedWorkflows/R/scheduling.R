# needs msLevel, precursorMZ, retentionTime
#' @export
header_to_schedule <- function(data) {
  data %>%
    filter(msLevel == 2) %>%
    group_by(precursorMZ, add = TRUE) %>%
    summarize(
      start = min(retentionTime),
      end = max(retentionTime)
    ) %>%
    ungroup(precursorMZ)
}

# needs ms_level, pre_mz, rt
#' @export
header_to_schedule_clean <- function(data) {
  data %>%
    filter(ms_level == 2) %>%
    group_by(pre_mz, add = TRUE) %>%
    summarize(
      start = min(rt),
      end = max(rt)
    ) %>%
    ungroup(pre_mz)
}

# needs start, end
#' @export
schedule_to_count <- function(data) {
  data %>%
    tidyr::pivot_longer(c(start, end), names_to = "event", values_to = "rt") %>%
    arrange(rt) %>%
    mutate(count = cumsum(c(start = 1, end = -1)[event]))
}

#' @export
integrate_over_count <- function(rt, count) {
  sum((tail(rt, -1) - head(rt, -1)) * head(count, -1))
}

# needs windows
#' @export
gen_optimal_split <- function(windows, windows_irt, split_count, iterations = 10000) {
  best_int <- Inf
  best_split <- NULL

  windows_irt_split <-
    windows_irt %>%
    mutate(split = list(seq_len(split_count))) %>%
    tidyr::unnest(split)

  for (i in 1:iterations) {
    curr_split <- sample(rep(1:split_count, length.out = nrow(windows)))


    split_df <-
      windows %>%
      group_by(split = curr_split) %>%
      tidyr::unnest(windows) %>%
      bind_rows(windows_irt_split) %>%
      schedule_to_count() %>%
      summarize(int = integrate_over_count(rt, count ^ 2))

    curr_sum <- sum(split_df$int)

    if (curr_sum < best_int) {
      best_int <- curr_sum
      best_split <- curr_split
    }
  }

  list(
    split = best_split,
    windows =
      windows %>%
      group_by(split = best_split) %>%
      tidyr::unnest(windows) %>%
      bind_rows(windows_irt_split)
  )
}

