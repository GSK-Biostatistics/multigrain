# .graph_optimise_ga returns valid object

    Code
      ga_res <- .graph_optimise_ga(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, nsim = ctrl$nsim_global,
        global_opts = ctrl$global_opt)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

---

    Code
      ga_res <- .graph_optimise_ga(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, global_opts = ctrl$
          global_opt, num_threads = cran_cores(), nsim = ctrl$nsim_global)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

# .graph_optimise_local returns valid object

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, local_opts = ctrl$
          local_opt)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

---

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, local_opts = ctrl$
          local_opt, num_threads = cran_cores())
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

# x0 passed successfully from .graph_optimise_ga to local

    Code
      ga_x0 <- as.vector(.graph_optimise_ga(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, global_opts = ctrl$
          global_opt, nsim = ctrl$nsim_global)$ga_output@solution[1, ])
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

---

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, local_opts = ctrl$
          local_opt, x0 = ga_x0)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

---

    Code
      ga_x0 <- as.vector(.graph_optimise_ga(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, global_opts = ctrl$
          global_opt, num_threads = 1L, nsim = ctrl$nsim_global)$ga_output@solution[1,
        ])
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

---

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, alpha = alpha, local_opts = ctrl$
          local_opt, num_threads = 1L, x0 = ga_x0)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

# Optimise 4m conjunctive power: local search only

    Code
      graph_optimise_4m_result <- graph_optimise(pvals = pvals_4m, alpha = alpha,
        graph_constraint = no_constr, trial_success = conjunctive_4m_power, control = ctrl,
        global_search = FALSE)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# Optimise 4m conjunctive power: include global search

    Code
      graph_optimise_4m_result <- graph_optimise(pvals = pvals_4m, alpha = alpha,
        graph_constraint = no_constr, trial_success = conjunctive_4m_power, control = ctrl)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# Optimise 4m conjunctive power: include global search & verbose

    Code
      graph_optimise(pvals = pvals_4m, alpha = alpha, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, control = ctrl, verbose = "foo")
    Condition
      Error in `graph_optimise()`:
      ! `verbose` must be a logical vector, not the string "foo".

# Optimise 2m: E2E check

    Code
      graph_2m <- graph_optimize(pvals_2, graph_constraint = gc, trial_success = ts,
        alpha = 0.025, control = ctrl)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# graph_optimise with num_threads: local search only

    Code
      graph_optimise_4m_result <- graph_optimise(pvals = pvals[, 1:4], alpha = alpha,
      graph_constraint = no_constr, trial_success = conjunctive_4m_power,
      num_threads = cran_cores(), control = ctrl, global_search = FALSE)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# graph_optimise with num_threads: include global search

    Code
      graph_optimise_4m_result <- graph_optimise(pvals = pvals[, 1:4], alpha = alpha,
      graph_constraint = no_constr, trial_success = conjunctive_4m_power,
      num_threads = cran_cores(), control = ctrl)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# graph_optimise parallel vs serial identical results (1 thread)

    Code
      result_serial <- graph_optimise(pvals = pvals_3m, alpha = alpha,
        graph_constraint = no_constr, trial_success = disjunctive_3m_power, control = ctrl,
        global_search = FALSE)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

---

    Code
      result_parallel_1t <- graph_optimise(pvals = pvals_3m, alpha = alpha,
        graph_constraint = no_constr, trial_success = disjunctive_3m_power,
        num_threads = 1L, control = ctrl, global_search = FALSE)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# graph_optimise: 2m E2E check with parallelization

    Code
      graph_2m <- graph_optimise(pvals_2m, graph_constraint = gc, trial_success = ts,
        alpha = 0.025, num_threads = 2L, control = ctrl)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# graph_optimise stores settings correctly

    Code
      result <- graph_optimise(pvals = pvals_3m, alpha = alpha, graph_constraint = no_constr,
        trial_success = disjunctive_3m_power, num_threads = cran_cores(), control = ctrl,
        global_search = FALSE)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# graph_optimise $power reflects pruned graph, not pre-pruned

    Code
      result <- graph_optimise(pvals = pvals[, 1:6], alpha = 0.025, graph_constraint = graph_constraint_free(
        6), trial_success = avg_6m_power, control = ctrl)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      

# deprecation message for power_nsim_local

    Code
      graph_optimise_4m_result <- graph_optimise(pvals = pvals, alpha = alpha,
        graph_constraint = no_constr, trial_success = conjunctive_4m_power,
        power_nsim_local = 200)
    Condition
      Error:
      ! `power_nsim_local` was deprecated in multigrain 0.2.0 and is now defunct.
      i Please use `control_nsim_local()` instead.

# deprecation message for power_nsim_global

    Code
      graph_optimise_4m_result <- graph_optimise(pvals = pvals, alpha = alpha,
        graph_constraint = no_constr, trial_success = conjunctive_4m_power,
        power_nsim_global = 200)
    Condition
      Error:
      ! `power_nsim_global` was deprecated in multigrain 0.2.0 and is now defunct.
      i Please use `control_nsim_global()` instead.

