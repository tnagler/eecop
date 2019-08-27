get_psi_expectile <- function(y) {
    function(theta, t = 0.5)
        t * (y - theta) * (y >= theta) - (1 - t) * (theta - y) * (y < theta)
}

get_psi_quantile <- function(y) {
    function(theta, t = 0.5)
        t * (y >= theta) - (1 - t) * (y < theta)
}
