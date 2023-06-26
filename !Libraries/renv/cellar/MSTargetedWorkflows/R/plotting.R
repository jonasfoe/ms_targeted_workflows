#' @export
normalize_neg_scale <- function(plot, yvals, mode = c("both", "top", "bottom"), breaks_per_axis = 3, top_scale_factor = 1) {
  mode = rlang::arg_match(mode)

  max_yval <- max(yvals, 0)
  if (max_yval == 0)
    max_yval <- 1

  min_yval <- min(yvals, 0)
  if (min_yval == 0)
    min_yval <- -1

  if (mode == "top")
    min_yval_vis <- 0
  else
    min_yval_vis <- min_yval

  if (mode == "bottom")
    max_yval_vis <- 0
  else
    max_yval_vis <- max_yval

  fact <- max_yval / -min_yval / top_scale_factor

  tr <- function(x) {
    pos <- tidyr::replace_na(x < 0, FALSE)
    x[pos] <- x[pos] * fact
    x
  }
  tr_i <- function(x) {
    pos <- tidyr::replace_na(x < 0, FALSE)
    x[pos] <- x[pos] / fact
    x
  }

  breaks_upper <- scales::extended_breaks(breaks_per_axis)(c(0, max_yval))
  breaks_upper <- breaks_upper[breaks_upper != 0]
  breaks_lower <- scales::extended_breaks(breaks_per_axis)(c(0, min_yval))
  breaks_lower <- breaks_lower[breaks_lower != 0]
  breaks <- c(breaks_lower, 0, breaks_upper)
  labeller <- function(x) scales::scientific(abs(x))

  plot <- plot + scale_y_continuous(trans = scales::trans_new(name = "negdiff", transform = tr, inverse = tr_i), breaks = breaks, labels = labeller)

  plot$coordinates$limits$y <- c(min_yval_vis, max_yval_vis)
  # plot <- suppressWarnings(plot + coord_cartesian(xlim = plot$coordinates$limits$x, ylim = c(min_yval_vis, max_yval_vis)))

  plot
}

#' @export
plot_chrom <- function(data, rt_lims, rt_shift_ref = 0, show_points = FALSE, ions, legend_position = c("inside", "outside"), ...) {
  legend_position <- rlang::arg_match(legend_position)

  if (nrow(data) > 0)
    assert_that(all(data$i >= 0))

  data_plot <- data

  #   if (!missing(rt_lims))
  #     data_plot <- filter(data_plot, rt > rt_lims[1], rt < rt_lims[2])
  if (!missing(ions))
    data_plot <- filter(data_plot, transition_name %in% ions)

  assert_that(nrow(data_plot) > 0)

  data_plot <- mutate(data_plot, is_pre = transition == "precursor")

  if (has_name(data_plot, "transition_rank_lib"))
    data_plot <- arrange(data_plot, !is_pre, transition_rank_lib)
  else
    data_plot <- arrange(data_plot, !is_pre, desc(transition_ordinal), transition_mz)

  data_plot <- mutate(data_plot, transition_name = factor(transition_name, levels = unique(transition_name)))

  plot <-
    data_plot %>%
    ggplot(aes(rt, i, color = transition_name)) +
    geom_line()

  if (show_points)
    plot <- plot + geom_point()

  plot <-
    plot +
    geom_hline(yintercept = 0) +
    scale_y_continuous(labels = scales::scientific) +
    labs(x = "Retention Time (min)", y = "Intensity (a.u.)", color = "Ion")

  if (!missing(rt_lims))
    plot <- plot + coord_cartesian(xlim = rt_lims)

  if (legend_position == "inside")
    plot <- plot + theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.title.align = 1
    )

  plot
}

#' @export
plot_chrom_vs_ref <- function(data_upper, data_lower, rt_lims, i_lims, rt_shift_ref = 0, show_points = FALSE, ions, legend_position = c("inside", "outside"), ...) {
  legend_position <- rlang::arg_match(legend_position)

  if (nrow(data_upper) > 0)
    assert_that(all(data_upper$i >= 0))

  if (nrow(data_lower) > 0) {
    data_lower <- mutate(data_lower, i = -i, rt = rt + rt_shift_ref)
    assert_that(all(data_lower$i <= 0))
  }

  if (nrow(data_upper) > 0 && nrow(data_lower) == 0)
    mode <- "top"
  else if (nrow(data_upper) == 0 && nrow(data_lower) > 0)
    mode <- "bottom"
  else
    mode <- "both"

  data_plot <- bind_rows(sample = data_upper, reference = data_lower, .id = "type")

  #   if (!missing(rt_lims))
  #     data_plot <- filter(data_plot, rt > rt_lims[1], rt < rt_lims[2])
  if (!missing(ions))
    data_plot <- filter(data_plot, transition_name %in% ions)

  assert_that(nrow(data_plot) > 0)

  data_plot <- mutate(data_plot, is_pre = transition == "precursor")

  if (has_name(data_plot, "transition_rank_lib"))
    data_plot <- arrange(data_plot, !is_pre, transition_rank_lib)
  else
    data_plot <- arrange(data_plot, !is_pre, desc(transition_ordinal), transition_mz)

  data_plot <- mutate(data_plot, transition_name = factor(transition_name, levels = unique(transition_name)))

  plot <-
    data_plot %>%
    tidyr::unite("group", type, transition_name, transition_mz, remove = FALSE) %>%
    ggplot(aes(rt, i, group = group, color = transition_name)) +
    geom_line()

  if (show_points)
    plot <- plot + geom_point()

  plot <-
    plot +
    geom_hline(yintercept = 0) +
    labs(x = "Retention Time (min)", y = "Intensity (a.u.)", color = "Ion")


  if (!missing(rt_lims))
    plot <- plot + coord_cartesian(xlim = rt_lims)

  plot_i_lims <- range(data_plot$i)
  if (!missing(i_lims)) {
    if (!is.na(i_lims[1])) plot_i_lims[1] <- i_lims[1]
    if (!is.na(i_lims[2])) plot_i_lims[2] <- i_lims[2]
  }

  plot <- normalize_neg_scale(plot, yvals = plot_i_lims, mode = mode, ...)

  if (legend_position == "inside")
    plot <- plot + theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.title.align = 1
    )

  plot
}

#' @export
plot_chrom_vs_ref_by_heavy <- function(data, heavy_col, ...) {
  plot_chrom_vs_ref(
    data %>% filter(!(!!enquo(heavy_col))),
    data %>% filter(!!enquo(heavy_col)),
    ...
  )
}

#' @export
plot_chrom_vs_ref_by_col <- function(data, upper_col, ...) {
  plot_chrom_vs_ref(
    data %>% filter(!!enquo(upper_col)),
    data %>% filter(!(!!enquo(upper_col))),
    ...
  )
}

#' @export
plot_chrom_dual <- function(data, upper_col, rt_lims, rt_shift_ref = 0, show_points = FALSE, legend_position = c("outside", "inside"), ...) {
  legend_position <- rlang::arg_match(legend_position)
  upper_col <- enquo(upper_col)
  data_upper <- data %>% filter(!!upper_col)
  data_lower <- data %>% filter(!(!!upper_col))

  if (nrow(data_upper) > 0)
    assert_that(all(data_upper$i >= 0))

  if (nrow(data_lower) > 0) {
    data_lower <- mutate(data_lower, i = -i, rt = rt + rt_shift_ref)
    assert_that(all(data_lower$i <= 0))
  }

  if (nrow(data_upper) > 0 && nrow(data_lower) == 0)
    mode <- "top"
  else if (nrow(data_upper) == 0 && nrow(data_lower) > 0)
    mode <- "bottom"
  else
    mode <- "both"

  data_plot <- bind_rows(sample = data_upper, reference = data_lower, .id = "type")

  data_plot <- mutate(data_plot, is_pre = transition == "precursor")

  if (has_name(data_plot, "transition_rank_lib"))
    data_plot <- arrange(data_plot, !is_pre, transition_rank_lib, desc(transition_ordinal), transition_mz)
  else
    data_plot <- arrange(data_plot, !is_pre, desc(transition_ordinal), transition_mz)

  data_plot <- mutate(data_plot, transition_name = factor(transition_name, levels = unique(transition_name)))

  plot <-
    data_plot %>%
    tidyr::unite("group", type, transition_name, transition_mz, remove = FALSE) %>%
    ggplot(aes(rt, i, group = group, color = transition_name)) +
    geom_line()

  if (show_points)
    plot <- plot + geom_point()

  plot <-
    plot +
    geom_hline(yintercept = 0) +
    labs(x = "Retention Time (min)", y = "Intensity (a.u.)", color = "Ion")


  if (!missing(rt_lims))
    plot <- plot + coord_cartesian(xlim = rt_lims)

  plot <- normalize_neg_scale(plot, yvals = data_plot$i, mode = mode, ...)

  if (legend_position == "inside")
    plot <- plot + theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.title.align = 1
    )

  plot
}

#' @export
plot_peptide_sequence_coverage <- function(sequence, pre_charge, transition_types, transition_ordinals, underscore_positions = integer(), letter_spacing = 0.8, mod_fontsize_factor = 0.8, mod_color = "grey30", width = 7, height = 0.65, fontsize = 6.5, vjust = 0, units = c("in", "cm", "mm")) {
  units <- rlang::arg_match(units)
  assert_that(
    is.string(sequence), is.count(pre_charge), is.number(letter_spacing), is.number(width), is.number(height), is.number(fontsize), is.number(vjust),
    is.character(transition_types), rlang::is_integerish(transition_ordinals), all(transition_ordinals > 0), length(transition_types) == length(transition_ordinals),
    rlang::is_integerish(underscore_positions), all(underscore_positions > 0)
  )

  unit_factors <- c(`in` = 1, `cm` = 2.54, `mm` = 25.4)
  unit_factor <- unit_factors[units]

  # 1 unit on the x-axis will correspond to the width of 'A'
  letter_a_width <- graphics::strwidth("A", family = "sans", units = "inches")
  seq_position_text <- first(str_match_all(sequence, "([A-Z])(\\[[^\\]]*\\])*"))

  # figure out the placement of the sequence text
  text_df <-
    tibble(
      position_text = na.omit(c(t(seq_position_text[, 2:3]))),
      is_pos = str_detect(position_text, "^[A-Z]$"),
      i_pos = replace(cumsum(is_pos), !is_pos, NA_integer_),
      x_width = graphics::strwidth(position_text, family = "sans", units = "inches") * ifelse(is_pos, 1, mod_fontsize_factor) / letter_a_width,
      x_spacing_post = c((diff(is_pos) >= 0), FALSE) * letter_spacing,
      x_end = cumsum(x_width + x_spacing_post) - x_spacing_post / 2,
      x_pos = x_end - x_spacing_post / 2 - x_width / 2,
      x_start = lag(x_end, n = 1, default = 0)
    )
  # relocate(.after = last_col(), x_start, x_pos, x_end)

  seq_length <- sum(text_df$is_pos)

  transitions_df <-
    tibble(transition_type = transition_types, transition_ordinal = transition_ordinals) %>%
    mutate(
      # which series does the transition belong to?
      transition_y_pos = case_when(
        transition_type == "precursor" ~ 0,
        transition_type %in% c("x", "y", "z") ~ 1,
        transition_type %in% c("a", "b", "c") ~ -1
      ),
      # the x location where the fragmentation cut is placed
      transition_x_pre_seq_pos = case_when(
        transition_type == "precursor" ~ 1,
        transition_type %in% c("x", "y", "z") ~ seq_length - transition_ordinal + 1,
        transition_type %in% c("a", "b", "c") ~ transition_ordinal + 1
      )
    ) %>%
    # the count annotation for the fragmentation
    count(transition_x_pre_seq_pos, transition_y_pos) %>%
    # the x placement of the cut
    left_join(select(text_df, i_pos, transition_x_pos = x_start), by = join_by(transition_x_pre_seq_pos == i_pos))

  shapes_template_df <- tribble(
    ~type, ~i, ~x, ~y,
    "precursor", 1, 0.25, -0.45,
    "precursor", 2, -0.1, -0.45,
    "precursor", 3, -0.1, 0.45,
    "precursor", 4, 0.25, 0.45,
    "lower", 1, -0.8, -0.5,
    "lower", 2, 0, -0.5,
    "lower", 3, 0, 0.1,
    "upper", 1, 0, -0.1,
    "upper", 2, 0, 0.5,
    "upper", 3, 0.8, 0.5,
    "both", 1, -0.8, -0.5,
    "both", 2, 0, -0.5,
    "both", 3, 0, 0.5,
    "both", 3, 0.8, 0.5,
  )

  shapes_df <-
    transitions_df %>%
    reframe(.by = transition_x_pos, type = case_when(
      all(transition_y_pos == 0) ~ "precursor",
      n() == 2 && !any(transition_y_pos == 0) ~ "both",
      n() == 1 && transition_y_pos > 0 ~ "upper",
      n() == 1 && transition_y_pos < 0 ~ "lower"
    )) %>%
    tibble::rowid_to_column("shape_i") %>%
    left_join(shapes_template_df, by = join_by(type), relationship = "many-to-many")

  assert_that(!anyNA(shapes_df$type))

  underscores_df <-
    text_df %>%
    filter(i_pos %in% underscore_positions) %>%
    select(position_text, i_pos, x_width, x_pos)

  # the range of the axes that enables consistent placements of all elements when applied via coord_cartesian
  xspan <- 38.35 / unit_factor * width / fontsize
  yspan <- 20 / unit_factor * height / fontsize
  # By design, all the content fits into this y limit (neg. to pos.)
  ylimits_content <- 0.95
  # if a lot of space is available, a lot of the yspan remains blank
  yspan_blank <- yspan - ylimits_content * 2

  geom_cutpath <- partial(geom_path, color = "red", linewidth = fontsize * .07)
  geom_cuttext <- partial(geom_text, color = "red", size = fontsize * .4)

  text_df %>%
    ggplot(aes(x = x_pos)) +
    geom_cutpath(data = shapes_df, aes(x = transition_x_pos + x, y = y, group = shape_i)) +
    geom_segment(data = underscores_df, aes(x = x_pos - x_width / 2, xend = x_pos + x_width / 2), y = -0.4, yend = -0.4, linewidth = fontsize * .07) +
    # sequence text
    geom_text(
      aes(
        label = position_text,
        # the letter baseline appears to be a ca. -0.283, which is where the mod text would end up at size 0
        y = ifelse(is_pos, 0, -0.283 * (1 - !!mod_fontsize_factor)),
        size = ifelse(is_pos, !!fontsize, !!fontsize * !!mod_fontsize_factor),
        color = ifelse(is_pos, "black", !!mod_color)
      ),
      vjust = .5, family = "sans"
    ) +
    # charge state +'s
    annotate("text", label = str_dup("+", pre_charge), x = last(text_df$x_end), y = .35,
      vjust = .5, hjust = 0, family = "sans", size = fontsize * .6
    ) +
    # counter for transitions
    geom_cuttext(data = dplyr::filter(transitions_df, transition_y_pos != 0),
      aes(label = n, x = transition_x_pos + .4 * transition_y_pos, y = .6 * transition_y_pos, vjust = transition_y_pos / -2 + .5),
      hjust = .5
    ) +
    # counter for precursor
    geom_cuttext(data = dplyr::filter(transitions_df, transition_y_pos == 0),
      aes(label = n, x = transition_x_pos - .2, y = transition_y_pos), vjust = .5, hjust = 1
    ) +
    scale_size_identity() +
    scale_color_identity() +
    coord_cartesian(
      xlim = c(-.8, xspan - .8),
      # ylim = c(-ylimits_content - yspan_blank / 2, ylimits_content + yspan_blank / 2,
      ylim = c(-ylimits_content - (yspan_blank) * (1 - vjust), ylimits_content + (yspan_blank) * vjust),
      expand = FALSE
    ) +
    theme_void()
}
