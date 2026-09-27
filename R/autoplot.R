# ggplot2 versions of the package's plots.
#
# ggplot2 is optional (Suggests). The autoplot() methods are registered
# with ggplot2's generic when the package is loaded (see .onLoad in
# zzz.R); theme_dsge() and the colour scales fail with an informative
# error when ggplot2 is not installed. Each method builds a tidy data
# frame and returns a ggplot object, so users can add layers, change
# facets or apply their own theme.

.dsge_need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This function needs the ggplot2 package: ",
         "install.packages(\"ggplot2\")", call. = FALSE)
  }
  invisible(TRUE)
}

#' ggplot2 Theme and Colour Scales for dsge Plots
#'
#' `theme_dsge()` is the ggplot2 theme used by the package's [autoplot()]
#' methods: a light horizontal grid, no panel borders, bold left-aligned
#' titles and facet labels, and the legend on top.
#' `scale_colour_dsge()` and `scale_fill_dsge()` apply the package's
#' categorical palette, whose order is checked to stay distinguishable
#' under the common forms of colour vision deficiency.
#'
#' These functions need the \pkg{ggplot2} package.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param ... Passed to [ggplot2::scale_colour_manual()] or
#'   [ggplot2::scale_fill_manual()].
#'
#' @return A ggplot2 theme or scale.
#'
#' @seealso [autoplot-dsge] for the ggplot2 versions of the package's
#'   plots.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   df <- data.frame(t = rep(1:20, 2), y = c(0.8^(0:19), 0.5^(0:19)),
#'                    s = rep(c("a", "b"), each = 20))
#'   ggplot2::ggplot(df, ggplot2::aes(t, y, colour = s)) +
#'     ggplot2::geom_line() +
#'     scale_colour_dsge() +
#'     theme_dsge()
#' }
#'
#' @export
theme_dsge <- function(base_size = 11, base_family = "") {
  .dsge_need_ggplot()
  el_text <- ggplot2::element_text
  el_line <- ggplot2::element_line
  blank <- ggplot2::element_blank()
  ggplot2::theme_minimal(base_size = base_size,
                         base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = blank,
      panel.grid.major.x = blank,
      panel.grid.major.y = el_line(colour = .DSGE_INK_GRID,
                                   linewidth = 0.4),
      axis.ticks.x = el_line(colour = .DSGE_INK_TICK, linewidth = 0.4),
      axis.text = el_text(colour = .DSGE_INK_AXIS),
      axis.title = el_text(colour = .DSGE_INK_AXIS,
                           size = ggplot2::rel(0.9)),
      strip.text = el_text(face = "bold", hjust = 0,
                           colour = .DSGE_INK_TITLE,
                           size = ggplot2::rel(1.0)),
      plot.title = el_text(face = "bold", hjust = 0,
                           colour = .DSGE_INK_TITLE),
      plot.subtitle = el_text(hjust = 0, colour = .DSGE_INK_AXIS),
      plot.title.position = "plot",
      legend.position = "top",
      legend.title = blank,
      legend.text = el_text(colour = .DSGE_INK_AXIS),
      panel.spacing = ggplot2::unit(1.1, "lines")
    )
}

#' @rdname theme_dsge
#' @export
scale_colour_dsge <- function(...) {
  .dsge_need_ggplot()
  ggplot2::scale_colour_manual(values = .dsge_palette(16L), ...)
}

#' @rdname theme_dsge
#' @export
scale_color_dsge <- scale_colour_dsge

#' @rdname theme_dsge
#' @export
scale_fill_dsge <- function(...) {
  .dsge_need_ggplot()
  ggplot2::scale_fill_manual(values = .dsge_palette(16L), ...)
}


#' ggplot2 Versions of the dsge Plots
#'
#' When \pkg{ggplot2} is installed, [ggplot2::autoplot()] draws the
#' package's results as ggplot objects, which can be extended with further
#' layers, re-faceted or restyled like any other ggplot. They mirror the
#' base-graphics [plot()] methods and use [theme_dsge()].
#'
#' @param object The result to plot: impulse responses ([irf()],
#'   [irf_2nd_order()]), a forecast ([forecast()]), a variance
#'   decomposition ([variance_decomposition()]), a historical shock
#'   decomposition ([shock_decomposition()]), smoothed states
#'   ([smooth_states()]), a perfect-foresight path ([perfect_foresight()]),
#'   an OccBin simulation ([simulate_occbin()]) or a Bayesian fit
#'   ([bayes_dsge()]).
#' @param impulse,response Shocks and responses to show (default: all, less
#'   auxiliary lag states of imported Dynare models and shocks that enter
#'   as variables).
#' @param ci Show confidence bands when available.
#' @param variables Variables to show.
#' @param which States or observables to show (names or indices).
#' @param type For a Bayesian fit, `"trace"` or `"density"`.
#' @param pars Parameters to show for a Bayesian fit.
#' @param ... Unused.
#'
#' @return A ggplot object.
#'
#' @name autoplot-dsge
#' @aliases autoplot.dsge_irf autoplot.dsge_irf_2nd autoplot.dsge_forecast
#'   autoplot.dsge_variance_decomposition autoplot.dsge_decomposition
#'   autoplot.dsge_smoothed autoplot.dsge_perfect_foresight
#'   autoplot.dsge_occbin autoplot.dsge_bayes
#' @usage
#' \method{autoplot}{dsge_irf}(object, impulse = NULL, response = NULL,
#'   ci = TRUE, ...)
#' \method{autoplot}{dsge_irf_2nd}(object, variables = NULL, ...)
#' \method{autoplot}{dsge_forecast}(object, ...)
#' \method{autoplot}{dsge_variance_decomposition}(object, ...)
#' \method{autoplot}{dsge_decomposition}(object, which = NULL, ...)
#' \method{autoplot}{dsge_smoothed}(object, which = NULL, ...)
#' \method{autoplot}{dsge_perfect_foresight}(object, variables = NULL, ...)
#' \method{autoplot}{dsge_occbin}(object, variables = NULL, ...)
#' \method{autoplot}{dsge_bayes}(object, type = c("trace", "density"),
#'   pars = NULL, ...)
#'
#' @seealso [theme_dsge()], [plot.dsge_irf()]
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   m <- dsge_model(
#'     obs(p   ~ beta * lead(p) + kappa * x),
#'     unobs(x ~ lead(x) - (r - lead(p) - g)),
#'     obs(r   ~ psi * p + u),
#'     state(u ~ rhou * u),
#'     state(g ~ rhog * g),
#'     fixed = list(beta = 0.99),
#'     start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
#'   )
#'   sol <- solve_dsge(m, params = c(kappa = 0.1, psi = 1.5, rhou = 0.7,
#'                                   rhog = 0.9))
#'   ggplot2::autoplot(irf(sol, periods = 16))
#'   ggplot2::autoplot(variance_decomposition(sol, horizon = c(1, 4, 8)))
#' }
NULL

autoplot.dsge_irf <- function(object, impulse = NULL, response = NULL,
                              ci = TRUE, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  dat <- object$data
  if (!is.null(impulse)) dat <- dat[dat$impulse %in% impulse, ]
  if (!is.null(response)) {
    dat <- dat[dat$response %in% response, ]
  } else if (nrow(dat) > 0L) {
    dat <- dat[dat$response %in% .dsge_irf_default_responses(dat), ]
  }
  if (nrow(dat) == 0L) {
    stop("Nothing to plot: no matching impulse/response pairs.",
         call. = FALSE)
  }
  # drop responses that are zero for every shock
  scale <- max(abs(dat$value), na.rm = TRUE)
  nonzero <- tapply(abs(dat$value), dat$response, max, na.rm = TRUE) >
    1e-10 * max(scale, 1e-300)
  keep <- names(nonzero)[nonzero]
  if (length(keep) > 0L) dat <- dat[dat$response %in% keep, ]

  dat$response <- factor(dat$response, levels = unique(dat$response))
  dat$impulse <- factor(dat$impulse, levels = unique(dat$impulse),
                        labels = paste("Shock:", unique(dat$impulse)))
  has_ci <- ci && all(c("lower", "upper") %in% names(dat)) &&
    !all(is.na(dat$lower))
  n_imp <- nlevels(dat$impulse)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$period,
                                         y = .data$value)) +
    ggplot2::geom_hline(yintercept = 0, colour = .DSGE_INK_REF,
                        linewidth = 0.4)
  if (has_ci) {
    lvl <- if (!is.null(object$level)) object$level else 0.95
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper,
                                        fill = sprintf("%g%% band",
                                                       100 * lvl)),
                           alpha = 0.25) +
      ggplot2::scale_fill_manual(values = .DSGE_INK_PRIMARY)
  }
  p <- p +
    ggplot2::geom_line(colour = .DSGE_INK_PRIMARY, linewidth = 0.9) +
    ggplot2::labs(x = "Periods after shock", y = NULL) +
    theme_dsge()
  if (n_imp > 1L) {
    # one row per shock, each panel on its own scale; strips read
    # "shock -> response"
    strip <- function(labels) {
      list(paste0(sub("^Shock: ", "", labels$impulse), " \u2192 ",
                  labels$response))
    }
    p + ggplot2::facet_wrap(ggplot2::vars(.data$impulse, .data$response),
                            scales = "free_y",
                            ncol = nlevels(dat$response),
                            labeller = strip)
  } else {
    p + ggplot2::facet_wrap(ggplot2::vars(.data$response),
                            scales = "free_y") +
      ggplot2::labs(title = paste("Responses to a shock to",
                                  sub("^Shock: ", "", levels(dat$impulse))))
  }
}

autoplot.dsge_irf_2nd <- function(object, variables = NULL, ...) {
  ir <- structure(list(
    data = data.frame(period = object$period, impulse = object$shock,
                      response = object$variable, value = object$response,
                      stringsAsFactors = FALSE)),
    class = "dsge_irf")
  autoplot.dsge_irf(ir, response = variables, ci = FALSE)
}

autoplot.dsge_forecast <- function(object, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  fc <- object$forecasts
  vars <- unique(fc$variable)
  has_sd <- "sd" %in% names(fc) && !all(is.na(fc$sd))
  has_hist <- !is.null(object$history) && nrow(object$history) > 0
  n_hist <- if (has_hist) nrow(object$history) else 0L
  hist_show <- if (has_hist) max(1L, n_hist - 3L * object$horizon + 1L)
               else 1L

  fc$t <- n_hist + fc$period
  lines <- data.frame(t = fc$t, variable = fc$variable, value = fc$value,
                      series = "Forecast", stringsAsFactors = FALSE)
  if (has_hist) {
    idx <- hist_show:n_hist
    hist <- do.call(rbind, lapply(vars, function(v) {
      data.frame(t = idx, variable = v, value = object$history[idx, v],
                 series = "Data", stringsAsFactors = FALSE)
    }))
    # join the forecast to the last observation
    last <- hist[hist$t == n_hist, ]
    last$series <- "Forecast"
    lines <- rbind(hist, last, lines)
  }
  lines$series <- factor(lines$series, levels = c("Data", "Forecast"))

  p <- ggplot2::ggplot()
  if (has_hist) {
    p <- p + ggplot2::annotate("rect", xmin = n_hist + 0.5, xmax = Inf,
                               ymin = -Inf, ymax = Inf, fill = "#F4F3EF")
  }
  if (has_sd) {
    z <- c("95%" = stats::qnorm(0.975), "80%" = stats::qnorm(0.9),
           "50%" = stats::qnorm(0.75))
    bands <- do.call(rbind, lapply(names(z), function(l) {
      b <- data.frame(t = fc$t, variable = fc$variable,
                      lower = fc$value - z[[l]] * fc$sd,
                      upper = fc$value + z[[l]] * fc$sd,
                      band = l, stringsAsFactors = FALSE)
      if (has_hist) {
        last <- data.frame(t = n_hist, variable = vars,
                           lower = object$history[n_hist, vars],
                           upper = object$history[n_hist, vars],
                           band = l, stringsAsFactors = FALSE)
        b <- rbind(last, b)
      }
      b
    }))
    bands$band <- factor(bands$band, levels = c("50%", "80%", "95%"))
    p <- p +
      ggplot2::geom_ribbon(data = bands,
                           ggplot2::aes(x = .data$t, ymin = .data$lower,
                                        ymax = .data$upper,
                                        fill = .data$band)) +
      ggplot2::scale_fill_manual(values = c(
        "50%" = paste0(.DSGE_INK_PRIMARY, "66"),
        "80%" = paste0(.DSGE_INK_PRIMARY, "40"),
        "95%" = paste0(.DSGE_INK_PRIMARY, "24")))
  }
  p +
    ggplot2::geom_line(data = lines,
                       ggplot2::aes(x = .data$t, y = .data$value,
                                    colour = .data$series,
                                    linewidth = .data$series)) +
    ggplot2::scale_colour_manual(values = c(Data = .DSGE_INK_TITLE,
                                            Forecast = .DSGE_INK_PRIMARY)) +
    ggplot2::scale_linewidth_manual(values = c(Data = 0.5, Forecast = 0.9),
                                    guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y",
                        ncol = if (length(vars) > 3L) 2L else 1L) +
    ggplot2::labs(x = "Period", y = NULL) +
    theme_dsge()
}

autoplot.dsge_variance_decomposition <- function(object, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  pct_labels <- function(x) paste0(x, "%")
  if (object$type == "unconditional") {
    pct <- object$contribution_pct
    dat <- data.frame(
      variable = factor(rep(rownames(pct), ncol(pct)),
                        levels = rev(rownames(pct))),
      shock = factor(rep(colnames(pct), each = nrow(pct)),
                     levels = colnames(pct)),
      share = as.vector(pct))
    return(
      ggplot2::ggplot(dat, ggplot2::aes(x = .data$share,
                                        y = .data$variable,
                                        fill = .data$shock)) +
        ggplot2::geom_col(width = 0.7, colour = "white",
                          linewidth = 0.3,
                          position = ggplot2::position_stack(reverse = TRUE)) +
        ggplot2::scale_x_continuous(labels = pct_labels,
                                    expand = c(0, 0)) +
        scale_fill_dsge() +
        ggplot2::labs(title = "Unconditional variance decomposition",
                      x = "Share of variance", y = NULL) +
        theme_dsge() +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_line(
                         colour = .DSGE_INK_GRID, linewidth = 0.4))
    )
  }
  arr <- object$contribution_pct          # n_h x n_o x n_e
  dat <- expand.grid(h = seq_along(object$horizon),
                     o = seq_along(object$obs_names),
                     e = seq_along(object$shock_names))
  dat$share <- arr[cbind(dat$h, dat$o, dat$e)]
  dat$horizon <- factor(object$horizon[dat$h], levels = object$horizon)
  dat$variable <- factor(object$obs_names[dat$o],
                         levels = object$obs_names)
  dat$shock <- factor(object$shock_names[dat$e],
                      levels = object$shock_names)
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$horizon, y = .data$share,
                                    fill = .data$shock)) +
    ggplot2::geom_col(width = 0.75, colour = "white", linewidth = 0.3,
                      position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_y_continuous(labels = pct_labels, expand = c(0, 0)) +
    scale_fill_dsge() +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable)) +
    ggplot2::labs(title = "Forecast-error variance decomposition",
                  x = "Horizon (periods)", y = NULL) +
    theme_dsge()
}

autoplot.dsge_decomposition <- function(object, which = NULL, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  dec <- object$decomposition
  n_T <- dim(dec)[1]
  n_comp <- dim(dec)[3]
  obs <- object$obs_names
  if (is.null(which)) which <- seq_along(obs)
  if (is.character(which)) which <- match(which, obs)
  comps <- c(object$shock_names[seq_len(n_comp - 1L)], "Initial conditions")
  dat <- expand.grid(t = seq_len(n_T), o = which, k = seq_len(n_comp))
  dat$value <- dec[cbind(dat$t, dat$o, dat$k)]
  dat$variable <- factor(obs[dat$o], levels = obs[which])
  dat$component <- factor(comps[dat$k], levels = comps)
  total <- stats::aggregate(value ~ t + variable, data = dat, FUN = sum)
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$t, y = .data$value)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$component), width = 0.85,
                      position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_hline(yintercept = 0, colour = .DSGE_INK_REF,
                        linewidth = 0.4) +
    ggplot2::geom_line(data = total, colour = .DSGE_INK_TITLE,
                       linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = c(.dsge_palette(n_comp - 1L),
                                          .DSGE_MUTED)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), ncol = 1,
                        scales = "free_y") +
    ggplot2::labs(title = "Historical shock decomposition",
                  subtitle = "line: the series; bars: contribution of each shock",
                  x = "Period", y = NULL) +
    theme_dsge()
}

autoplot.dsge_smoothed <- function(object, which = NULL, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  st <- object$smoothed_states
  nm <- object$state_names
  if (is.null(which)) which <- seq_len(min(ncol(st), 9L))
  if (is.character(which)) which <- match(which, nm)
  dat <- do.call(rbind, lapply(which, function(i) {
    d <- data.frame(t = seq_len(nrow(st)), state = nm[i], value = st[, i],
                    stringsAsFactors = FALSE)
    if (!is.null(object$smoothed_states_var)) {
      s <- sqrt(pmax(object$smoothed_states_var[, i], 0))
      d$lower <- d$value - 2 * s
      d$upper <- d$value + 2 * s
    }
    d
  }))
  dat$state <- factor(dat$state, levels = nm[which])
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$t, y = .data$value)) +
    ggplot2::geom_hline(yintercept = 0, colour = .DSGE_INK_REF,
                        linewidth = 0.4)
  if ("lower" %in% names(dat)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = .DSGE_INK_PRIMARY, alpha = 0.2)
  }
  p +
    ggplot2::geom_line(colour = .DSGE_INK_PRIMARY, linewidth = 0.7) +
    ggplot2::facet_wrap(ggplot2::vars(.data$state), scales = "free_y") +
    ggplot2::labs(title = "Smoothed states",
                  subtitle = if ("lower" %in% names(dat))
                    "deviation from steady state, ± 2 s.d." else
                    "deviation from steady state",
                  x = "Period", y = NULL) +
    theme_dsge()
}

.dsge_paths_long <- function(mat, vars, label) {
  data.frame(t = rep(seq_len(nrow(mat)), length(vars)),
             variable = rep(vars, each = nrow(mat)),
             value = as.vector(mat[, vars, drop = FALSE]),
             series = label, stringsAsFactors = FALSE)
}

autoplot.dsge_perfect_foresight <- function(object, variables = NULL, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  all_data <- cbind(object$controls, object$states)
  if (is.null(variables)) {
    active <- apply(abs(all_data), 2, max) > 1e-10
    variables <- .dsge_drop_aux(colnames(all_data)[active],
                                object$shock_names)
  }
  dat <- .dsge_paths_long(all_data, variables, "path")
  dat$variable <- factor(dat$variable, levels = variables)
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$t, y = .data$value)) +
    ggplot2::geom_hline(yintercept = 0, colour = .DSGE_INK_REF,
                        linewidth = 0.4) +
    ggplot2::geom_line(colour = .DSGE_INK_PRIMARY, linewidth = 0.9) +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y") +
    ggplot2::labs(title = "Perfect-foresight path",
                  subtitle = "deviation from steady state",
                  x = "Period", y = NULL) +
    theme_dsge()
}

autoplot.dsge_occbin <- function(object, variables = NULL, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  con <- cbind(object$controls, object$states)
  unc <- cbind(object$controls_unc, object$states_unc)
  con_vars <- vapply(object$constraints, function(c) c$variable,
                     character(1))
  if (is.null(variables)) {
    active <- colnames(con)[apply(abs(con), 2, max) > 1e-10]
    variables <- unique(c(con_vars, active))
  }
  dat <- rbind(.dsge_paths_long(con, variables, "With constraint"),
               .dsge_paths_long(unc, variables, "Without constraint"))
  dat$variable <- factor(dat$variable, levels = variables)
  bind <- do.call(rbind, lapply(seq_along(con_vars), function(ci) {
    b <- which(object$binding[, ci])
    if (length(b) == 0L || !con_vars[ci] %in% variables) return(NULL)
    data.frame(variable = factor(con_vars[ci], levels = variables),
               xmin = b - 0.5, xmax = b + 0.5)
  }))
  p <- ggplot2::ggplot()
  if (!is.null(bind)) {
    p <- p + ggplot2::geom_rect(
      data = bind, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
      ymin = -Inf, ymax = Inf, fill = .DSGE_INK_SECONDARY, alpha = 0.15)
  }
  p +
    ggplot2::geom_hline(yintercept = 0, colour = .DSGE_INK_REF,
                        linewidth = 0.4) +
    ggplot2::geom_line(data = dat,
                       ggplot2::aes(x = .data$t, y = .data$value,
                                    colour = .data$series,
                                    linetype = .data$series),
                       linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = c(
      "With constraint" = .DSGE_INK_PRIMARY,
      "Without constraint" = .DSGE_INK_NEUTRAL)) +
    ggplot2::scale_linetype_manual(values = c(
      "With constraint" = "solid", "Without constraint" = "dashed")) +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y") +
    ggplot2::labs(title = "Occasionally binding constraint",
                  subtitle = if (!is.null(bind))
                    "shaded: periods in which the constraint binds",
                  x = "Period", y = NULL) +
    theme_dsge()
}

autoplot.dsge_bayes <- function(object, type = c("trace", "density"),
                                pars = NULL, ...) {
  .dsge_need_ggplot()
  .data <- ggplot2::.data
  type <- match.arg(type)
  post <- object$posterior                 # iter x par x chain
  par_names <- dimnames(post)[[2]]
  if (!is.null(pars)) {
    bad <- setdiff(pars, par_names)
    if (length(bad) > 0L) {
      stop("Unknown parameter(s): ", paste(bad, collapse = ", "),
           call. = FALSE)
    }
  } else {
    pars <- par_names
  }
  n_it <- dim(post)[1]
  n_ch <- dim(post)[3]
  dat <- expand.grid(iter = seq_len(n_it), p = match(pars, par_names),
                     chain = seq_len(n_ch))
  dat$value <- post[cbind(dat$iter, dat$p, dat$chain)]
  dat$parameter <- factor(par_names[dat$p], levels = pars)
  dat$chain <- factor(paste("Chain", dat$chain))

  if (type == "trace") {
    return(
      ggplot2::ggplot(dat, ggplot2::aes(x = .data$iter, y = .data$value,
                                        colour = .data$chain)) +
        ggplot2::geom_line(linewidth = 0.3, alpha = 0.8) +
        scale_colour_dsge() +
        ggplot2::facet_wrap(ggplot2::vars(.data$parameter), ncol = 1,
                            scales = "free_y") +
        ggplot2::labs(title = "MCMC traces", x = "Iteration", y = NULL) +
        theme_dsge()
    )
  }
  # prior densities over each parameter's posterior range
  priors <- object$priors
  prior_df <- do.call(rbind, lapply(pars, function(pn) {
    j <- match(pn, par_names)
    if (is.null(priors) || length(priors) < j || is.null(priors[[j]])) {
      return(NULL)
    }
    r <- range(post[, j, ])
    xs <- seq(r[1], r[2], length.out = 200)
    dens <- vapply(xs, function(v) exp(dprior(priors[[j]], v)), numeric(1))
    if (!any(is.finite(dens) & dens > 0)) return(NULL)
    data.frame(parameter = pn, value = xs, density = dens)
  }))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density(ggplot2::aes(linetype = "Posterior"),
                          fill = .DSGE_INK_PRIMARY, alpha = 0.2,
                          colour = .DSGE_INK_PRIMARY, linewidth = 0.7)
  if (!is.null(prior_df)) {
    prior_df$parameter <- factor(prior_df$parameter, levels = pars)
    p <- p + ggplot2::geom_line(data = prior_df,
                                ggplot2::aes(y = .data$density,
                                             linetype = "Prior"),
                                colour = .DSGE_INK_NEUTRAL,
                                linewidth = 0.6)
  }
  p +
    ggplot2::scale_linetype_manual(values = c(Posterior = "solid",
                                              Prior = "dashed")) +
    ggplot2::facet_wrap(ggplot2::vars(.data$parameter), scales = "free") +
    ggplot2::labs(title = "Posterior densities", x = NULL, y = NULL) +
    theme_dsge()
}
