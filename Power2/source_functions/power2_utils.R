
#  ---------------------------------------
# barplot of cluster prevalences
barplot_clusters <- function(cluster_prev) {
  cluster_prev %>%
    select(cluster, p) %>%
    mutate(q = 1 - p) %>%
    pivot_longer(-cluster) %>%
    ggplot() + theme_bw() +
    geom_bar(aes(x = as.factor(cluster), y = (1 - value) * 100, fill = name), position = "stack", stat = "identity") +
    scale_fill_manual(values = c(grey(0.8), "coral1"), guide = "none") +
    xlab("Cluster") + ylab("Prevalance")
}

#  ---------------------------------------
# draw from beta-binomial model with a defined level of intra-cluster
# correlation
draw_overdispersed <- function(n_clusters, n_samp, prev, ICC) {
  s1 <- prev * (1/ICC - 1)
  s2 <- (1 - prev) * (1/ICC - 1)
  cluster_p <- rbeta(n_clusters, shape1 = s1, shape2 = s2)
  ret <- rbinom(n_clusters, n_samp, cluster_p)
  return(ret)
}

#  ---------------------------------------
# draw from over-dispersed model and calculate CIs by simple-random-sampling,
# and taking into account design effect. Plot both CIs for comparison.
compare_CIs <- function(n_clusters, n_samp, prev, ICC) {
  
  p <- draw_overdispersed(n_clusters, n_samp, prev, ICC) / n_samp
  p_bar <- mean(p)
  var_clust <- var(p) / (n_clusters - 1)
  var_srs <- p_bar*(1 - p_bar) / (n_clusters*n_samp - 1)
  Deff <- var_clust / var_srs
  SE <- sqrt(p_bar*(1 - p_bar) / (n_clusters*n_samp - 1))
  SE_adj <- sqrt(p_bar*(1 - p_bar) / (n_clusters*n_samp / Deff - 1))
  
  data.frame(type = c("Simple random\nsampling", "Cluster sampling"),
             p_bar = p_bar,
             lower = p_bar - 1.96*c(SE, SE_adj),
             upper = p_bar + 1.96*c(SE, SE_adj)) %>%
    ggplot() + theme_bw() +
    geom_hline(yintercept = 100*prev, linetype = "dashed") +
    geom_pointrange(aes(x = type, y = 100*p_bar, ymin = 100*lower, ymax = 100*upper)) +
    ylim(c(0, 100)) +
    xlab("") + ylab("Prevalence")
}

#  ---------------------------------------
# function to plot the distribution of the t-test statistic under the null and
# alternative hypotheses
plot_ttest_cluster <- function(p1, p2, n_samp, n_clusters, Deff, alpha = 0.05) {
  
  # check that p1 is largest
  if (p1 < p2) {
    stop("p1 must be larger than p2 by convention")
  }
  
  # calculate effective sample size
  n_eff <- n_samp / Deff
  
  # calculate the critical values under the null hypothesis
  critical_values <- qt(c(alpha / 2, 1 - alpha / 2), df = 2*n_eff - 2)
  
  # calculate the mean of the distribution under the alternative hypothesis
  alt_mean <- (p1 - p2) / sqrt( (p1*(1 - p1) + p2*(1 - p2)) / (n_eff*n_clusters) )
  
  # plot null and alternative distributions
  ret <- data.frame(x = seq(-10, 10, l = 201)) %>%
    mutate(t_dist_null = dt(x, df = 2*n_clusters - 2),
           t_dist_alt = dt(x, df = 2*n_clusters - 2, ncp = alt_mean)) %>%
    ggplot() + theme_bw() +
    geom_line(aes(x = x, y = t_dist_null)) +
    geom_line(aes(x = x, y = t_dist_alt), color = "red") +
    geom_vline(xintercept = critical_values, linetype = "dashed") +
    xlab("t statistic") + ylab("Probability")
  
  ret
}

#  ---------------------------------------
# returns the power under the cluster t-test
get_pow_ttest_cluster <- function(p1, p2, n_samp, n_clusters, Deff, alpha = 0.05) {
  n_eff <- n_samp / Deff
  alt_mean <- (p1 - p2) / sqrt( (p1*(1 - p1) + p2*(1 - p2)) / (n_eff*n_clusters) )
  pt(qt(alpha / 2, df = 2*n_clusters - 2), df = 2*n_clusters - 2, ncp = alt_mean) + 
    pt(qt(1 - alpha / 2, df = 2*n_clusters - 2), df = 2*n_clusters - 2, ncp = alt_mean, lower.tail = FALSE)
}

#  ---------------------------------------
# returns the maximum possible power under the cluster t-test (infinite sample
# size per cluster)
get_max_pow_ttest_cluster <- function(p1, p2, n_clusters, ICC, alpha = 0.05) {
  alt_mean <- (p1 - p2) / sqrt( (p1*(1 - p1) + p2*(1 - p2)) / (n_clusters / ICC) )
  pt(qt(alpha / 2, df = 2*n_clusters - 2), df = 2*n_clusters - 2, ncp = alt_mean) + 
    pt(qt(1 - alpha / 2, df = 2*n_clusters - 2), df = 2*n_clusters - 2, ncp = alt_mean, lower.tail = FALSE)
}

#  ---------------------------------------
# returns the power under the cluster one-sample t-test against a threshold
get_pow_ttest_thresh <- function(p, n_samp, n_clusters, ICC, p_thresh = 0.05, alpha = 0.05) {
  alt_mean <- (p - p_thresh) / sqrt( (1 + (n_samp - 1)*ICC) * p*(1 - p) / (n_samp*n_clusters) )
  pt(qt(alpha / 2, df = n_clusters - 1), df = n_clusters - 1, ncp = alt_mean) + 
    pt(qt(1 - alpha / 2, df = n_clusters - 1), df = n_clusters - 1, ncp = alt_mean, lower.tail = FALSE)
}

#  ---------------------------------------
# returns the maximum possible power under the one-sample cluster t-test
# (infinite sample size per cluster)
get_max_pow_ttest_thresh <- function(p, n_clusters, ICC, mu = 0.05, alpha = 0.05) {
  alt_mean <- (p - mu) / sqrt( ICC * p*(1 - p) / n_clusters )
  pt(qt(alpha / 2, df = n_clusters - 1), df = n_clusters - 1, ncp = alt_mean) + 
    pt(qt(1 - alpha / 2, df = n_clusters - 1), df = n_clusters - 1, ncp = alt_mean, lower.tail = FALSE)
}


