matrix_tolerance <- function(x, tolerance = 10e-13) {
    set.seed(1L)
    vec <- sample.int(100, size = x^2)

    temp_matrix <- matrix(
        vec,
        nrow = x,
        byrow = TRUE
    )

    diag(temp_matrix) <- 0

    sum_rows <- rowSums(temp_matrix)

    output_matrix <- (temp_matrix / sum_rows) + tolerance

    diag(output_matrix) <- 0

    output_matrix
}

matrix_tolerance_na <- function(x, num_nas = 2, tolerance = 10e-13) {
    set.seed(1L)
    vec <- sample.int(100, size = (x^2) - (num_nas * x))
    nas <- rep(NA, num_nas * x)

    mat_vec <- sample(c(vec, nas), size = x^2, replace = FALSE)

    temp_matrix <- matrix(
        mat_vec,
        nrow = x,
        byrow = TRUE
    )

    diag(temp_matrix) <- 0

    num_nas_per_row <- rowSums(is.na(temp_matrix))
    rows_with_one_na <- which(num_nas_per_row == 1)

    for (i in rows_with_one_na) {
        temp_matrix[i, ] <- tidyr::replace_na(temp_matrix[i, ], 0)
    }

    sum_rows <- rowSums(temp_matrix, na.rm = TRUE)

    output_matrix <- (temp_matrix / sum_rows) + tolerance

    diag(output_matrix) <- 0

    output_matrix
}
