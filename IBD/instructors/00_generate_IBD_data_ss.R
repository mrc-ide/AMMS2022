## .................................................................................
## Purpose: Generate data for AMMS IBD practical using `polySimIBD`
##
## Author: Nick Brazeau
##
## Date: 10 July, 2022
##
## Notes: Apologies, package is still undergoing changes and I may break it before AMMS.
## This commit works for our purposes
## .................................................................................
remotes::install_github("nickbrazeau/polySimIBD", ref = "develop")
library(polySimIBD)
library(vcfR)
library(tidyverse)
set.seed(48)

#............................................................
# Setup site demes based on practical figure
#...........................................................
# vectors must be ordered for population A, B, C, D, E
ds <- 12
demesizes <- rep(ds,2)
dwnsmpl <- split(cumsum(demesizes), f = 1:length(demesizes))
dwnsmpl <- sort(unlist(lapply(dwnsmpl, function(x,ds)sample((x-ds):x, size = 5), ds = ds)))
coimeans <- c(1,1)
m <- rep(0.25, 2)
# make symmetrical
migr_dist_mat <- matrix(c(1, 0,
                          2, 1),
                        ncol = 2, nrow = 2)
swfsim <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)),
                              N = demesizes,
                              m = m,
                              mean_coi = coimeans,
                              migr_dist_mat = migr_dist_mat,
                              rho = 1e-2,
                              tlim = 10)

#......................
# lift over into VCF
#......................
# get arg
ARG <- polySimIBD::get_arg(swfsim, host_index = dwnsmpl)
# convert to haploype matrix
hapmat <- polySimIBD::get_haplotype_matrix(ARG)
# Simulate Reads
this_coi <- swfsim$coi[dwnsmpl]
reads <- polySimIBD::sim_biallelic(COIs = this_coi,
                                   haplotypematrix = hapmat,
                                   shape1 = 1.544,
                                   shape2 = 0.620,
                                   coverage = 100,
                                   alpha = 1,
                                   overdispersion = 0.1,
                                   epsilon = 0.01)
# convert to VCF
make_vcf_from_reads <- function(reads, swfsim) {
  # checks
  # TODO

  # get FIXED
  fix <- tibble::tibble(
    CHROM = "CHROM1",
    POS = swfsim$pos,
    ID = 1:length(swfsim$pos),
    REF = "N",
    ALT = "N",
    QUAL = "100",
    FILTER = "PASS",
    INFO = "PolySimIBD=T"
  )

  # get gt
  assignGTfrombiWSNRAF <- function(wsnraf, cutoff = cutoff){
    if(wsnraf < cutoff) {
      return("0/0")
    } else if(wsnraf > (1-cutoff)) {
      return("1/1")
    } else {
      return("0/1")
    }
  }
  gtmat <- matrix(NA, nrow = nrow(reads$NRWSAcounts), ncol = ncol(reads$NRWSAcounts))
  for(i in 1:nrow(reads$NRWSAcounts)) {
    for (j in 1:ncol(reads$NRWSAcounts)) {
      gtmat[i,j] <- paste(c(
        assignGTfrombiWSNRAF(reads$NRWSAF[i,j], cutoff = 0.05),
        paste(c(reads$WS.coverage[i,j] - reads$NRWSAcounts[i,j], reads$NRWSAcounts[i,j]), collapse = ","),
        reads$WS.coverage[i,j]
      ), collapse = ":")
    }
  }


  # append format column and sample names
  gtmat <- cbind(FORMAT = "GT:AD:DP", gtmat)
  colnames(gtmat)[2:ncol(gtmat)] <- paste0("smpl", 1:ncol(reads$NRWSAcounts))

  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = "##fileformat=VCFv4.3", fix = as.matrix(fix), gt = gtmat)
  return(newvcfR)
}

simVCF <- make_vcf_from_reads(reads, swfsim)
# manually write over smpl names for participant convenience
colnames(simVCF@gt)[2:ncol(simVCF@gt)] <- c(paste0("A", 1:5), paste0("B", 1:5))
# save out
vcfR::write.vcf(simVCF, file = "data/simulated_sink_source_ibd.vcf.gz")


#............................................................
# playground
#...........................................................
# get MLE IBD
mipvcf <- MIPanalyzer::vcf2mipanalyzer_biallelic(vcfR = simVCF)
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


# historgram
mainplot <- ibd_long %>%
  ggplot() +
  geom_histogram(aes(x=malecotf, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "#d9d9d9") +
  xlab("IBD") + ylab("frequency (%)") +
  theme_classic()

insetplot <- ibd_long %>%
  ggplot() +
  geom_histogram(aes(x=malecotf, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "#d9d9d9") +
  xlab("IBD") + ylab("frequency (%)") +
  theme_classic() +
  coord_cartesian(xlim = c(0.5,1), ylim = c(0,0.8)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))
# cowplot
cowplot::ggdraw() +
  cowplot::draw_plot(mainplot, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(insetplot, x = 0.5, y= 0.3, width = 0.4, height = 0.4)


