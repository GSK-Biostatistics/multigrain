#' Extract power from rejection matrix.
#'
#' Reference implementation used for testing. Mirrors the output structure of
#' `calc_power_pvals()`.
#'
#' @param x  n × m logical/0-1 matrix; rows = simulations, cols = hypotheses.
#' @param f  (optional) list of custom functions, each taking a single vector of
#'           rejections and returning a scalar.
#' @returns A list with at least four elements:
#'   * `local_power`: column means
#'   * `exp_rejections`: expected # of rejections
#'   * `disj_power`: disjunctive power
#'   * `conj_power`: conjunctive power
#'
#'   Plus one element per function in `f`.
extract_power <- function(x, f = list()) {
    stopifnot(is.matrix(x), mode(x) %in% c("logical", "numeric"))

    m <- ncol(x)
    rej_counts <- rowSums(x)

    out <- list(
        local_power = colMeans(x),
        exp_rejections = mean(rej_counts),
        disj_power = mean(rej_counts > 0L),
        conj_power = mean(rej_counts == m)
    )

    if (is.function(f)) {
        f <- list(f)
    }

    if (length(f)) {
        nm <- names(f)
        if (is.null(nm)) {
            nm <- rep("", length(f))
        }
        empty <- nm == "" | is.na(nm)
        nm[empty] <- sprintf("func%d", which(empty))
        names(f) <- nm

        for (nm_i in nm) {
            out[[nm_i]] <- suppressWarnings(
                mean(
                    vapply(
                        seq_len(nrow(x)),
                        \(i) f[[nm_i]](x[i, ]),
                        FUN.VALUE = numeric(1L)
                    )
                )
            )
        }
    }

    out
}
