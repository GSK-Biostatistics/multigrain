# Truncated Cauchy perturbation centered on parent value.
# Draws from Cauchy(centre, scale) via rejection sampling.
.cauchy_perturb <- function(centre, lower, upper, scale = 1.0) {
    for (i in seq_len(100L)) {
        proposal <- centre + scale * tan(pi * (stats::runif(1) - 0.5))
        if (proposal >= lower && proposal <= upper) {
            return(proposal)
        }
    }
    # Fallback: uniform on [lower, upper]
    stats::runif(1, lower, upper)
}


# Factory: returns a multi-parameter Cauchy perturbation mutation closure.
#
# Perturbs each parameter within graph with probability p_param_mutate
#
# Designed for use with p_mutation = 1.0 (every individual mutated every
# generation), with exploitation handled by Nelder-Mead local search
# (optim = TRUE, poptim = 0.2) and elitism preserving the best solution(s).
#
# @param p_param_mutate Per-parameter mutation probability. Default 0.1.
# @param scale Cauchy scale parameter. Default 1.0 (matching ESCH).
# @returns A function(object, parent) for the [GA::ga()] mutation slot.
.make_cauchy_mutation_multi <- function(p_param_mutate = 0.1, scale = 1.0) {
    force(p_param_mutate)
    force(scale)
    function(object, parent) {
        parent_vec <- as.numeric(object@population[parent, ])
        d <- length(parent_vec)
        lower <- object@lower
        upper <- object@upper

        # Independently decide which parameters to mutate
        mutate_mask <- stats::runif(d) < p_param_mutate

        # Guarantee at least one mutation (avoid wasting a fitness evaluation)
        if (!any(mutate_mask)) {
            mutate_mask[sample.int(d, 1L)] <- TRUE
        }

        # Perturb selected parameters around their current values
        for (j in which(mutate_mask)) {
            parent_vec[j] <- .cauchy_perturb(
                parent_vec[j],
                lower[j],
                upper[j],
                scale
            )
        }

        parent_vec
    }
}
