#' @title shhh
#' @description sinks annoying calls to console
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' @title View recombination plot generato
#' @description stick breaking idea of recombination breakdown btwn
#' two samples that were originally clonal
view_recombo <- function(generations_apart) {
  # draw recombo
  matbk <- sort( runif(n = generations_apart, min = 0, max = 1) )
  patbk <- sort( runif(n = generations_apart, min = 0, max = 1) )
  # tidy and plot
  plot <- tibble::tibble(start = c(0,matbk, 0, patbk),
                         end = c(matbk,1, patbk, 1),
                         par = c(rep("mat", generations_apart+1),
                                 rep("pat", generations_apart+1)
                         )) %>%
    dplyr::mutate(
      ibd = ifelse(generations_apart != 0, rbinom(dplyr::n(), 1, 0.5), 0),
      ibdchar = ifelse(ibd == 0, "Orig", "New"),
      ymin = ifelse(par == "mat", 0.1, 1.1),
      ymax = ifelse(par == "mat", 1, 2)) %>%
    ggplot() +
    geom_rect(aes(xmin = start, xmax = end,
                  ymin = ymin, ymax = ymax,
                  fill = ibdchar)) +
    scale_fill_viridis_d("Identity") +
    labs(title = "Breakdown of Clonal Material over Generations",
         x = "Genomic Position") +
    theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
          axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
          legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "#000000", size = 1),
          axis.line.y = element_blank())

  # out
  return(plot)
}
