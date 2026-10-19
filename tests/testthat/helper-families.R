# Literal Algorithm 5.7 of Maurer and Bretz (2014), Bonferroni across families.
#
# Rejects one elementary hypothesis at a time and updates the family graph
# only when a family is emptied (all its members rejected). It is independent
# of the C++ kernel: adjusted p-values come from stats::p.adjust(), and the
# rejection order can be the lowest index or random. Returns an N x m logical
# rejection matrix. Also sourced by dev/review/families/verify.R.
alg57_oracle <- function(
    pvals,
    members,
    test,
    w,
    G,
    alpha = 0.025,
    order = c("lowest", "random")
) {
    order <- match.arg(order)
    test <- rep_len(test, length(members))
    family_of <- family_index0(members) + 1L

    rej <- matrix(FALSE, nrow = nrow(pvals), ncol = ncol(pvals))

    for (trial in seq_len(nrow(pvals))) {
        padj <- numeric(ncol(pvals))
        for (f in seq_along(members)) {
            idx <- members[[f]]
            padj[idx] <- stats::p.adjust(pvals[trial, idx], method = test[f])
        }
        rej[trial, ] <- alg57_trial(
            padj,
            members,
            family_of,
            w * alpha,
            G,
            order
        )
    }

    rej
}

# One trial of Algorithm 5.7 on within-family adjusted p-values `padj`, with
# initial local levels `a` and family transition matrix `g`.
alg57_trial <- function(padj, members, family_of, a, g, order) {
    active <- rep(TRUE, length(members))
    rejected <- rep(FALSE, length(padj))

    repeat {
        candidates <- which(!rejected & padj < a[family_of])
        if (length(candidates) == 0L) {
            break
        }
        j <- if (order == "lowest" || length(candidates) == 1L) {
            candidates[1]
        } else {
            sample(candidates, 1)
        }
        rejected[j] <- TRUE

        f <- family_of[j]
        if (all(rejected[members[[f]]])) {
            # Family f is emptied and leaves the graph
            active[f] <- FALSE
            updated <- bretz_update(a, g, f, others = which(active))
            a <- updated$a
            g <- updated$g
        }
    }

    rejected
}

# Bretz et al. (2009) update of local levels `a` and transition matrix `g`
# when node `f` leaves the graph; `others` are the nodes still in it.
bretz_update <- function(a, g, f, others) {
    a[others] <- a[others] + a[f] * g[f, others]

    g_new <- g
    for (i in others) {
        for (k in others) {
            denom <- 1 - g[i, f] * g[f, i]
            g_new[i, k] <- if (i == k || denom == 0) {
                0
            } else {
                (g[i, k] + g[i, f] * g[f, k]) / denom
            }
        }
    }
    g_new[f, ] <- 0
    g_new[, f] <- 0
    a[f] <- 0

    list(a = a, g = g_new)
}

# Random partition of 1:m into num_fam non-empty families; members are
# usually not contiguous columns.
random_partition <- function(m, num_fam) {
    fam <- sample(c(seq_len(num_fam), sample.int(num_fam, m - num_fam, TRUE)))
    unname(split(seq_len(m), factor(fam, levels = seq_len(num_fam))))
}

# Random family graph on num_fam nodes. With `full_rows = FALSE` each row of G
# is scaled to sum to less than 1.
random_family_graph <- function(num_fam, full_rows = TRUE) {
    if (num_fam == 1L) {
        return(list(w = 1, G = matrix(0, 1, 1)))
    }
    gr <- graph_random(num_fam)
    G <- unname(gr$trans_matrix)
    if (!full_rows) {
        G <- G * stats::runif(num_fam, 0.3, 0.9)
    }
    list(w = unname(gr$hyp_weight), G = G)
}

# 0-based family index of each column, as the kernels expect.
family_index0 <- function(members) {
    index <- integer(sum(lengths(members)))
    for (f in seq_along(members)) {
        index[members[[f]]] <- f - 1L
    }
    index
}

# Method codes used by adjust_pvals_families().
test_code <- function(test) {
    match(test, c("hommel", "hochberg", "holm")) - 1L
}

# Within-family adjustment through the C++ routine, taking R-style inputs.
adjust_families_r <- function(pvals, members, test) {
    adjust_pvals_families(
        pvals,
        family = family_index0(members),
        test = test_code(rep_len(test, length(members)))
    )
}

# Kernel on raw p-values, taking R-style inputs.
shortcut_families_r <- function(pvals, members, test, w, G, alpha = 0.025) {
    graph_shortcut_families(
        adjust_families_r(pvals, members, test),
        family = family_index0(members),
        alpha = alpha,
        w = w,
        G = G
    )
}
