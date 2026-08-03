# Determine if CRAN limits are active

cran_cores <- function() {
    check_cran <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    is_cran_check <- nzchar(check_cran) && check_cran == "TRUE"

    if (is_cran_check) {
        return(2L)
    }

    num_cores <- parallel::detectCores(logical = FALSE)

    if (is.na(num_cores)) {
        num_cores <- 1L
    }

    max(1L, min(7L, num_cores))
}
