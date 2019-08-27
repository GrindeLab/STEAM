get_thresh_analytic_new <- function(g, map, alpha = 0.05, type = "pval", searchint = c(1.96,8)){
  nchr <- length(unique(map$chr))
  delt.mat <- lapply(1:nchr, function(x) diff(map$cM[map$chr == x]))
  delt <- mean(unlist(delt.mat))
  L <- sum(unlist(delt.mat))
  Z <- uniroot(OU_approx, interval = searchint, beta = 0.01 * 
                 g, Delta = delt, length = L, chr = nchr, center = alpha, 
               test = "two-sided")$root
  if (type == "pval") {
    p <- 2 * pnorm(Z, lower.tail = F)
    return(p)
  }
  else if (type == "stat") {
    return(Z)
  }
  else {
    cat("Error: please specify type='stat' or 'pval'.")
  }
}