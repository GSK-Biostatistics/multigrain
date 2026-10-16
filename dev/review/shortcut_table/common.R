# Shared helpers for the shortcut table-kernel review scripts.
#
# Each script takes the library holding the build under test as its first
# argument, so v1.0.0 and the fast-shortcut build can be installed side by
# side (e.g. R CMD INSTALL --library=<lib>) and compared in separate R
# processes.

load_build <- function(lib) {
    suppressPackageStartupMessages(library(multigrain, lib.loc = lib))
    desc <- utils::packageDescription("multigrain", lib.loc = lib)
    build <- desc[["Config/multigrain/build"]]
    list(
        version = desc$Version,
        build = if (is.null(build)) desc$Version else build,
        lib = lib
    )
}

# Equicorrelated one-sided p-values (rho = 0.5). Means rise from 1.5 to 3.5
# across H1..H(m-1) and H(m) is null; "low" power shifts every non-null mean
# down by 1 and "very_high" shifts it up by 1.
sim_pvals <- function(N, m, power = c("high", "low", "very_high"), seed = 1) {
    power <- match.arg(power)
    set.seed(seed)
    means <- if (m == 1) 3 else c(seq(1.5, 3.5, length.out = m - 1), 0)
    if (power == "low") {
        means <- pmax(means - 1, 0)
    } else if (power == "very_high") {
        means[means > 0] <- means[means > 0] + 1
    }
    common <- stats::rnorm(N)
    z <- sqrt(0.5) * common + sqrt(0.5) * matrix(stats::rnorm(N * m), N, m)
    z <- sweep(z, 2, means, "+")
    stats::pnorm(z, lower.tail = FALSE)
}

# Dense random graph from graph_random() (rows of G sum to 1).
dense_graph <- function(m, seed) {
    set.seed(seed)
    gr <- multigrain::graph_random(m)
    list(w = unname(gr$hyp_weight), G = unname(gr$trans_matrix))
}

# Sparse graph: about half the off-diagonal edges and a third of the weights
# zeroed and renormalised (at least one survives); about one row in five is
# then scaled to sum to 0.7.
sparse_graph <- function(m, seed) {
    gr <- dense_graph(m, seed)
    set.seed(seed + 1e6)
    w <- gr$w
    if (m > 1) {
        drop <- stats::runif(m) < 1 / 3
        if (all(drop)) drop[sample.int(m, 1)] <- FALSE
        w[drop] <- 0
        w <- w / sum(w)
    }
    G <- gr$G
    for (i in seq_len(m)) {
        off <- setdiff(seq_len(m), i)
        if (length(off) > 1) {
            drop <- stats::runif(length(off)) < 0.5
            if (all(drop)) drop[sample.int(length(off), 1)] <- FALSE
            G[i, off[drop]] <- 0
        }
        if (sum(G[i, ]) > 0) G[i, ] <- G[i, ] / sum(G[i, ])
        if (stats::runif(1) < 0.2) G[i, ] <- G[i, ] * 0.7
    }
    list(w = w, G = G)
}

# Run graph_shortcut() with the kernel forced by option.
shortcut_forced <- function(kernel, pvals, alpha, w, G, num_threads = NULL) {
    old <- options(multigrain.shortcut_kernel = kernel)
    on.exit(options(old))
    if (is.null(num_threads)) {
        multigrain:::graph_shortcut(pvals, alpha, w, G)
    } else {
        multigrain:::graph_shortcut_parallel(
            pvals,
            alpha,
            w,
            G,
            num_threads = num_threads
        )
    }
}

# R replay of one trial through both kernels, recording every comparison
# p_j < c_j each kernel made. Used only to classify disagreements.
replay_trial <- function(p, w, G, alpha, crit_table) {
    m <- length(p)
    # Sequential kernel (first rejectable hypothesis, graph update per step)
    a <- w * alpha
    rej_seq <- rep(FALSE, m)
    cmp_seq <- NULL
    repeat {
        cmp_seq <- rbind(cmp_seq, cbind(j = which(!rej_seq), p = p[!rej_seq], c = a[!rej_seq]))
        i <- which(!rej_seq & p < a)[1]
        if (is.na(i)) break
        rej_seq[i] <- TRUE
        if (all(rej_seq)) break
        a <- a + a[i] * G[i, ]
        G_new <- matrix(0, m, m)
        for (r in seq_len(m)) {
            denom <- 1 - G[r, i] * G[i, r]
            if (denom != 0) {
                for (s in seq_len(m)) {
                    if (r != s) G_new[r, s] <- (G[r, s] + G[r, i] * G[i, s]) / denom
                }
            }
        }
        G <- G_new
        G[i, ] <- 0
        G[, i] <- 0
        a[i] <- 0
    }
    # Table kernel
    S <- rep(TRUE, m)
    cmp_tab <- NULL
    repeat {
        mask <- sum(2^(which(S) - 1))
        c_S <- crit_table[mask + 1, ]
        cmp_tab <- rbind(cmp_tab, cbind(j = which(S), p = p[S], c = c_S[S]))
        r <- S & p < c_S
        if (!any(r)) break
        S <- S & !r
    }
    gaps <- function(cmp) abs(cmp[, "p"] - cmp[, "c"]) / pmax(cmp[, "c"], .Machine$double.xmin)
    list(
        rej_seq = rej_seq,
        rej_tab = !S,
        min_rel_gap = min(c(gaps(cmp_seq), gaps(cmp_tab))),
        cmp_seq = cmp_seq,
        cmp_tab = cmp_tab
    )
}
