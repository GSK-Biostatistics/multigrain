graph_violation_score_r <- function(w, G) {
    w <- as.numeric(w)
    G <- as.matrix(G)
    v <- 0
    v <- v +
        sum(pmax(0, -w), na.rm = TRUE) +
        sum(pmax(0, w - 1), na.rm = TRUE)
    v <- v +
        sum(pmax(0, -G), na.rm = TRUE) +
        sum(pmax(0, G - 1), na.rm = TRUE)
    v <- v + abs(sum(w, na.rm = TRUE) - 1)
    v <- v + sum(abs(rowSums(G, na.rm = TRUE) - 1))
    v
}
