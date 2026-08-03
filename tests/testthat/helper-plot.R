equivalent_ggplot2 <- function(x, y) {
    tmp1 <- tempfile(fileext = ".svg")
    tmp2 <- tempfile(fileext = ".svg")

    suppressMessages(ggplot2::ggsave(tmp1, x))
    suppressMessages(ggplot2::ggsave(tmp2, y))

    tools::md5sum(tmp1) == tools::md5sum(tmp2)
}
