test_that("new_family_spec() can create an empty object", {
    empty_fs <- new_family_spec()

    expect_s3_class(empty_fs, "multigrain_family_spec")
    expect_named(empty_fs, c("members", "test"))
    expect_identical(family_spec_get_m(empty_fs), 0L)
})

test_that("family_spec() works", {
    fams <- family_spec(list(primary = 1, secondary = c(4, 2, 3)))

    expect_s3_class(fams, "multigrain_family_spec")
    expect_identical(
        fams$members,
        list(primary = 1L, secondary = 2:4)
    )
    expect_identical(fams$test, c(primary = "hommel", secondary = "hommel"))
    expect_identical(family_spec_get_m(fams), 4L)
})

test_that("family_spec() names unnamed families", {
    fams <- family_spec(list(1, 2:3))

    expect_named(fams$members, c("F1", "F2"))
    expect_named(fams$test, c("F1", "F2"))
})

test_that("family_spec() takes one test per family", {
    fams <- family_spec(list(1, 2:3, 4), test = c("holm", "hochberg", "hommel"))

    expect_identical(unname(fams$test), c("holm", "hochberg", "hommel"))
})

test_that("family_spec() validates `members`", {
    expect_snapshot(error = TRUE, {
        family_spec(1:3)
    })
    expect_snapshot(error = TRUE, {
        family_spec(list())
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(1, 2.5))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(1, integer(0), 2))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(0:1, 2))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(1:2, 2:3))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(1, 3))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(1, 4))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(a = 1, 2))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(a = 1, a = 2)) # nolint: duplicate_argument_linter.
    })
})

test_that("family_spec() validates `test`", {
    expect_snapshot(error = TRUE, {
        family_spec(list(1, 2), test = "bonferroni")
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(1, 2, 3), test = c("holm", "hommel"))
    })
    expect_snapshot(error = TRUE, {
        family_spec(list(1, 2), test = 1)
    })
})

test_that("family_spec_get_index() gives the family of each hypothesis", {
    fams <- family_spec(list(c(1, 4), 2, c(3, 5)))

    expect_identical(family_spec_get_index(fams), c(1L, 2L, 3L, 1L, 3L))
})

test_that("check_family_spec()", {
    fams <- family_spec(list(1, 2:3))

    expect_no_error(check_family_spec(fams))
    expect_no_error(check_family_spec(NULL, allow_null = TRUE))
    expect_snapshot(error = TRUE, {
        check_family_spec(list(1, 2:3))
    })
})

test_that("check_families_dims()", {
    fams <- family_spec(list(1, 2:3))
    pvals <- matrix(0.5, nrow = 2, ncol = 3)

    expect_no_error(check_families_dims(fams, pvals, num_nodes = 2))
    expect_snapshot(error = TRUE, {
        check_families_dims(fams, pvals[, 1:2], num_nodes = 2)
    })
    expect_snapshot(error = TRUE, {
        check_families_dims(fams, pvals, num_nodes = 3)
    })
})

test_that(".prepare_pvals() adjusts within families only when supplied", {
    pvals <- pvals_fixture[1:100, 1:4]
    fams <- family_spec(list(1, 2:4), test = c("hommel", "holm"))

    expect_identical(.prepare_pvals(pvals, NULL), pvals)

    expected <- pvals
    expected[, 2:4] <- t(apply(
        pvals[, 2:4],
        1,
        stats::p.adjust,
        method = "holm"
    ))
    expect_identical(.prepare_pvals(pvals, fams), unname(expected))
})

test_that("family_spec print and summary methods", {
    fams <- family_spec(
        list(primary = 1, secondary = 2:4),
        test = c("hommel", "hochberg")
    )

    expect_snapshot(print(fams))
    expect_snapshot(summary(fams))
})
