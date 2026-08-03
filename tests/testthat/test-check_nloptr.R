test_that("nloptr check functions", {
    # Rosenbrock Banana function and gradient in separate functions
    eval_f <- function(x) {
        100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
    }

    eval_grad_f <- function(x) {
        c(
            -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
            200 * (x[2] - x[1] * x[1])
        )
    }

    # initial values
    x0 <- c(-1.2, 1)

    opts <- list(
        algorithm = "NLOPT_LD_LBFGS",
        xtol_rel = 1.0e-8
    )

    # solve Rosenbrock Banana function
    res <- nloptr::nloptr(
        x0 = x0,
        eval_f = eval_f,
        eval_grad_f = eval_grad_f,
        opts = opts
    )

    expect_true(is_nloptr(res))
    expect_no_error(check_nloptr(res))

    expect_no_error(check_nloptr(NULL, allow_null = TRUE))

    expect_error(
        check_nloptr("foo"),
        '`"foo"` must be a nloptr object, not the string "foo"'
    )
})
