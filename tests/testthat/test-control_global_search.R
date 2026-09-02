test_that("control_global_search is deprecated", {
    ctrl <- multigrain_control()

    expect_snapshot(error = TRUE, {
        control_global_search(ctrl, FALSE)
    })
})
