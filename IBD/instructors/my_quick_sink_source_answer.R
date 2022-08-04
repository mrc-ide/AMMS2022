## .................................................................................
## Purpose: my quick answer to source sink challenge
##
## Author: Nick Brazeau
##
## Date: 31 July, 2022
##
## Notes:
## .................................................................................
library(tidyverse)
library(vcfR)
library(MIPanalyzer)

ssVCF <- read.vcfR("data/simulated_sink_source_ibd.vcf.gz")
# get MLE IBD
mipvcf <- MIPanalyzer::vcf2mipanalyzer_biallelic(vcfR = ssVCF)
ibd <- MIPanalyzer::inbreeding_mle(x = mipvcf,
                                   f = seq(0.01, 0.99, 0.01),
                                   ignore_het = FALSE)
diag(ibd$mle) <- 1
colnames(ibd$mle) <- rownames(ibd$mle) <- colnames(simVCF@gt)[2:ncol(simVCF@gt)]

ibd_long <- broom::tidy(as.dist(t(ibd$mle))) %>%  # note, returns upper triangle
  magrittr::set_colnames(c("p1", "p2", "malecotf"))


# see if connections are as expected
library(tidygraph)
library(ggraph)
adj_graph <- ibd_long %>%
  tidygraph::as_tbl_graph(., directed = F)

ibd_long %>%
  tidygraph::as_tbl_graph(., directed = F) %>%
  dplyr::mutate(community = as.factor(tidygraph::group_louvain(weights = malecotf))) %>%
  tidygraph::activate("edges") %>%
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = malecotf,
                             color = malecotf)) +
  ggraph::geom_node_point(aes(color = community),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  ggraph::geom_node_text(aes(label = name), repel = T) +
  scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma") +
  scale_color_brewer(palette = "Set1") +
  ggraph::theme_graph() +
  theme(legend.position = "bottom")




ibd_long %>%
  dplyr::filter( (stringr::str_detect(p1, "A") & stringr::str_detect(p2, "A")) |
                   (stringr::str_detect(p1, "B") & stringr::str_detect(p2, "B"))
  ) %>%
  dplyr::mutate(
    pop = ifelse(stringr::str_detect(p1, "A"), "A", "B")) %>%
  dplyr::group_by(pop) %>%
  dplyr::summarise(
    withinIBD = mean(malecotf)
  )


