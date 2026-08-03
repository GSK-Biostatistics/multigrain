# .graph_optimise_ga returns valid object

    Code
      ga_res <- .graph_optimise_ga(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, nsim = ctrl$nsim_global, global_opts = ctrl$
          global_opt)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

---

    Code
      ga_res <- .graph_optimise_ga(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, global_opts = ctrl$global_opt,
        num_threads = cores, nsim = ctrl$nsim_global)
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

# .graph_optimise_local returns valid object

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, local_opts = ctrl$local_opt)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

---

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, local_opts = ctrl$local_opt,
        num_threads = cores)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

# x0 passed successfully from .graph_optimise_ga to local

    Code
      ga_x0 <- as.vector(.graph_optimise_ga(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, global_opts = ctrl$global_opt, nsim = ctrl$
          nsim_global)$ga_output@solution[1, ])
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

---

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_rand, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, local_opts = ctrl$local_opt, x0 = ga_x0)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

---

    Code
      ga_x0 <- as.vector(.graph_optimise_ga(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, global_opts = ctrl$global_opt,
        num_threads = 1L, nsim = ctrl$nsim_global)$ga_output@solution[1, ])
    Message
      i Running global optimization
      v Running global optimization [10ms]
      
      i Evaluating trial success of globally optimised graph
      v Evaluating trial success of globally optimised graph [10ms]
      

---

    Code
      loc_res <- .graph_optimise_local(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, local_opts = ctrl$local_opt,
        num_threads = 1L, x0 = ga_x0)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      

# Optimise 4m conjunctive power: local search only

    Code
      graph_optimise_4m_result <- graph_optimise(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, control = ctrl, global_search = FALSE)
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
      graph_optimise_4m_result <- graph_optimise(pvals = pvals_4m, graph_constraint = no_constr,
        trial_success = conjunctive_4m_power, control = ctrl)
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
      graph_optimise(pvals = pvals_4m, graph_constraint = no_constr, trial_success = conjunctive_4m_power,
        control = ctrl, verbose = "foo")
    Condition
      Error in `graph_optimise()`:
      ! `verbose` must be one of "info", "detail", or "silent", not "foo".

# Optimise 2m: E2E check

    Code
      graph_2m <- graph_optimize(pvals_2, graph_constraint = gc, trial_success = ts,
        control = ctrl)
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
      graph_optimise_4m_result <- graph_optimise(pvals = pvals[, 1:4],
      graph_constraint = no_constr, trial_success = conjunctive_4m_power,
      num_threads = cores, control = ctrl, global_search = FALSE)
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
      graph_optimise_4m_result <- graph_optimise(pvals = pvals[, 1:4],
      graph_constraint = no_constr, trial_success = conjunctive_4m_power,
      num_threads = cores, control = ctrl)
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
      result_serial <- graph_optimise(pvals = pvals_3m, graph_constraint = no_constr,
        trial_success = disjunctive_3m_power, control = ctrl, global_search = FALSE)
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
      result_parallel_1t <- graph_optimise(pvals = pvals_3m, graph_constraint = no_constr,
        trial_success = disjunctive_3m_power, num_threads = 1L, control = ctrl,
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
      

# graph_optimise: 2m E2E check with parallelization

    Code
      graph_2m <- graph_optimise(pvals_2m, graph_constraint = gc, trial_success = ts,
        num_threads = 2L, control = ctrl)
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
      result <- graph_optimise(pvals = pvals_3m, graph_constraint = no_constr,
        trial_success = disjunctive_3m_power, num_threads = cores, control = ctrl,
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
      result <- graph_optimise(pvals = pvals[, 1:6], graph_constraint = graph_constraint_free(
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
      

# graph_optimise complains when users pass ...

    Code
      graph_optim <- graph_optimise(pvals_4m, no_constr, conjunctive_4m_power, alpha = 0.025,
        control = ctrl)
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
      

---

    Code
      graph_optim <- graph_optimise(pvals_4m, no_constr, conjunctive_4m_power, alpha)
    Condition
      Error in `graph_optimise()`:
      ! `...` must be empty.
      x Problematic argument:
      * ..1 = alpha
      i Did you forget to name an argument?

---

    Code
      graph_optim <- graph_optimise(pvals_4m, no_constr, conjunctive_4m_power, ctrl)
    Condition
      Error in `graph_optimise()`:
      ! `...` must be empty.
      x Problematic argument:
      * ..1 = ctrl
      i Did you forget to name an argument?

# graph_optimise with old, logical verbose

    Code
      graph_optimise(pvals = pvals, graph_constraint = graph_constraint_free(4),
      trial_success = ts, num_threads = cores, global_search = FALSE, verbose = TRUE)
    Message
      i Running local optimization
      v Running local optimization [10ms]
      
      i Evaluating trial success of locally optimised graph
      v Evaluating trial success of locally optimised graph [10ms]
      
      i Pruning redundant weights and edges
      v Pruning redundant weights and edges [10ms]
      
      i Evaluating trial success of pruned graph
      v Evaluating trial success of pruned graph [10ms]
      
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
          H1     H2     H3     H4 
      0.2847 0.2204 0.2428 0.2521 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.5773 0.3182 0.1045
      H2 0.3585 0.0000 0.3021 0.3393
      H3 0.2456 0.3489 0.0000 0.4055
      H4 0.5453 0.1802 0.2745 0.0000
      
      Trial success function:
      r1 + r2 + r3 + r4
      
      Value of trial success measure:
      3.0458

---

    Code
      graph_optimise(pvals = pvals, graph_constraint = graph_constraint_free(4),
      trial_success = ts, num_threads = cores, global_search = FALSE, verbose = FALSE)
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
          H1     H2     H3     H4 
      0.2847 0.2204 0.2428 0.2521 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.5773 0.3182 0.1045
      H2 0.3585 0.0000 0.3021 0.3393
      H3 0.2456 0.3489 0.0000 0.4055
      H4 0.5453 0.1802 0.2745 0.0000
      
      Trial success function:
      r1 + r2 + r3 + r4
      
      Value of trial success measure:
      3.0458

