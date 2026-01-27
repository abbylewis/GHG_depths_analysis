pseudo_log_breaks <- function(n = 5, sigma = 1, base = 10) {
  force_all(n, base)
  n_default <- n
  function(x, n = n_default) {
    rng <- pseudo_log(range(x, na.rm = TRUE), sigma = sigma, base = base)
    min <- floor(rng[1])
    max <- ceiling(rng[2])

    if (max == min) {
      return(base^min)
    }

    by <- floor((max - min) / n) + 1
    breaks <- base^seq(min, max, by = by)
    relevant_breaks <- base^rng[1] <= breaks & breaks <= base^rng[2]
    if (sum(relevant_breaks) >= (n - 2)) {
      return(breaks)
    }

    # the easy solution to get more breaks is to decrease 'by'
    while (by > 1) {
      by <- by - 1
      breaks <- base^seq(min, max, by = by)
      relevant_breaks <- base^rng[1] <= breaks & breaks <= base^rng[2]
      if (sum(relevant_breaks) >= (n - 2)) {
        return(breaks)
      }
    }
    log_sub_breaks(rng, n = n, base = base)
  }
}
