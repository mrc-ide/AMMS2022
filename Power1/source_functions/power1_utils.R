
#  ---------------------------------------
# function to plot the distribution of the t-test statistic under the null and
# alternative hypotheses
plot_ttest <- function(d, s, n, alpha = 0.05) {
  
  # calculate the critical values under the null hypothesis
  critical_values <- qt(c(alpha / 2, 1 - alpha / 2), df = 2*n - 2)
  
  # calculate the mean of the distribution under the alternative hypothesis
  alt_mean <- d / sqrt( 2*s^2 / n)
  
  # plot null and alternative distributions
  ret <- data.frame(x = seq(-10, 10, l = 201)) %>%
    mutate(t_dist_null = dt(x, df = 2*n - 2),
           t_dist_alt = dt(x, df = 2*n - 2, ncp = alt_mean)) %>%
    ggplot() + theme_bw() +
    geom_line(aes(x = x, y = t_dist_null)) +
    geom_line(aes(x = x, y = t_dist_alt), color = "red") +
    geom_vline(xintercept = critical_values, linetype = "dashed") +
    xlab("t statistic") + ylab("Probability")
  
  ret
}

#  ---------------------------------------
# function to plot the distribution of the z-test statistic under the null and
# alternative hypotheses
plot_ztest <- function(p1, p2, n1, n2, alpha = 0.05) {
  
  # calculate the critical values under the null hypothesis
  critical_values <- qnorm(c(alpha / 2, 1 - alpha / 2))
  
  # calculate the mean of the distribution under the alternative hypothesis
  alt_mean <- (p1 - p2) / sqrt(p1*(1 - p1) / n1 + p2*(1 - p2) / n2)
  
  # plot null and alternative distributions
  ret <- data.frame(x = seq(-20, 20, l = 201)) %>%
    mutate(z_dist_null = dnorm(x),
           z_dist_alt = dnorm(x, mean = alt_mean)) %>%
    ggplot() + theme_bw() +
    geom_line(aes(x = x, y = z_dist_null)) +
    geom_line(aes(x = x, y = z_dist_alt), color = "red") +
    geom_vline(xintercept = critical_values, linetype = "dashed") +
    xlab("t statistic") + ylab("Probability")
  
  ret
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
