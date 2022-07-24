
popvar <- function(x) {
  mean((x - mean(x))^2)
}

sim_freqs <- function(N = 1e3,
                      demes = 1,
                      mig_rate = 0.01,
                      mut_rate = 1e-6,
                      loci = 1,
                      t_out = 0:(365 * 3),
                      plot_on = TRUE,
                      plot_homo = FALSE,
                      plot_expected_homo = FALSE) {
  
  # check inputs
  genescaper:::assert_single_pos_int(N, zero_allowed = FALSE)
  genescaper:::assert_single_pos_int(demes, zero_allowed = FALSE)
  genescaper:::assert_single_bounded(mig_rate)
  genescaper:::assert_single_bounded(mut_rate)
  genescaper:::assert_single_pos_int(loci, zero_allowed = FALSE)
  genescaper:::assert_vector_pos_int(t_out, zero_allowed = TRUE)
  genescaper:::assert_logical(plot_on)
  genescaper:::assert_logical(plot_homo)
  genescaper:::assert_logical(plot_expected_homo)
  
  # define migration matrix
  mig_mat <- matrix(mig_rate / demes, demes, demes)
  diag(mig_mat) <- 0
  diag(mig_mat) <- 1 - rowSums(mig_mat)
  
  # simulate from island model
  sim1 <- sim_wrightfisher(N = N,
                           L = loci,
                           alleles = 2,
                           mu = mut_rate,
                           t_out = t_out,
                           mig_mat = mig_mat,
                           initial_method = 2,
                           initial_params = rep(1e6, 2), silent = TRUE)
  
  # strip out unwated rows and columns, and calculate allele frequency
  ret <- sim1 %>%
    filter(allele == 1) %>%
    mutate(freq = count / N) %>%
    select(-allele, -count)
  
  # plots
  plot_ret <- NULL
  if (plot_on) {
    
    # plot frequencies
    plot1 <- ggplot(ret) + theme_bw() +
      geom_line(aes(x = time, y = freq, colour = interaction(locus, deme))) +
      ylim(c(0, 1)) +
      xlab("Time (days)") + ylab("Allele frequency") +
      scale_color_discrete(name = "Locus.Deme")
    
    # plot homozygosity
    if (plot_homo) {
      plot2 <- ggplot(ret) + theme_bw() +
        geom_line(aes(x = time, y = freq^2 + (1 - freq)^2, colour = interaction(locus, deme))) +
        ylim(c(0.5, 1)) +
        xlab("Time (days)") + ylab("Homozygosity") +
        scale_color_discrete(name = "Locus.Deme")
      
      # overlay expected homozygosity
      if (plot_expected_homo) {
        A <- 1/N + (1 - 1/N)*2*mut_rate
        B <- (1 - 4*mut_rate)*(1 - 1/N)
        J0 <- 0.5
        df_homo <- data.frame(x = t_out,
                              y = A*(B^(t_out + 1) - 1)/(B - 1) + (J0 - A)*B^t_out)
        
        plot2 <- plot2 +
          geom_line(aes(x = x, y = y), linetype = "dashed", data = df_homo)
      }
    }
    
    if (plot_homo) {
      plot_ret <- cowplot::plot_grid(plot1, plot2)
    } else {
      plot_ret <- plot1
    }
    print(plot_ret)
  }
  
  # return data and plot in list
  invisible(list(data = ret,
                 plot = plot_ret))
}
