
#  ---------------------------------------
# barplot of cluster prevalences
barplot_clusters <- function(cluster_prev) {
  cluster_prev %>%
    select(cluster, p) %>%
    mutate(q = 1 - p) %>%
    pivot_longer(-cluster) %>%
    ggplot() + theme_bw() +
    geom_bar(aes(x = cluster, y = (1 - value) * 100, fill = name), position = "stack", stat = "identity") +
    scale_fill_manual(values = c(grey(0.5), "firebrick2"), guide = "none") +
    xlab("Cluster") + ylab("Prevalance")
}

draw_overdispersed <- function(n_clusters, n_samp, prev, ICC) {
  s1 <- prev * (1/ICC - 1)
  s2 <- (1 - prev) * (1/ICC - 1)
  cluster_p <- rbeta(n_clusters, shape1 = s1, shape2 = s2)
  ret <- rbinom(n_clusters, n_samp, cluster_p)
  return(ret)
}

test1 <- function(n_clusters, n_samp, prev, ICC) {
  
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
# draw horizontal line using html
hrule <- function() {
  '<hr style="height:1px;border:none;color:#333;background-color:#333;">'
}

#  ---------------------------------------
# functions for beginning and ending an expandable answer
begin_button <- function(ID) {
  sprintf('<p><a class="btn-sm btn-primary" data-toggle="collapse" href="#collapseExample%s" role="button" aria-expanded="false" aria-controls="collapseExample%s">Click For Answer</a></p><div class="collapse" id="collapseExample%s"><div class="card card-body">', ID, ID, ID)
}
end_button <- function(ID) {
  '</div></div><br>'
}
