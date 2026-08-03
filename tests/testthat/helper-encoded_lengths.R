.expected_encoded_lengths <- function(constraint_w, constraint_G) {
    free_w <- sum(is.na(constraint_w))
    w_len <- max(free_w - 1L, 0L)
    row_free <- rowSums(is.na(constraint_G))
    g_len <- sum(pmax(row_free - 1L, 0L))

    list(
        free_w = as.integer(free_w),
        w_len = as.integer(w_len),
        g_len = as.integer(g_len)
    )
}
