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

#  ---------------------------------------
#' @title shhh
#' @description sinks annoying calls to console
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#  ---------------------------------------
#' @title COI Coin Toss Conceptual Visualization
#' @details for a given set of loci, plot results of binomial toss

COIn_toss <- function(COI = 1, loci = 3) {
  # assertions to do

  # realizations
  rt <- lapply(1:loci, function(x){rbinom(n = COI, size = 1, prob = 0.5)})
  # out
  dplyr::data_frame(rt) %>%
    dplyr::mutate(loci = paste0("Loci_", 1:dplyr::n())) %>%
    tidyr::unnest(., cols = "rt") %>%
    dplyr::mutate(coin = case_when(rt == 0 ~ "H",
                                   rt == 1 ~ "T")) %>%
    dplyr::group_by(loci) %>%
    dplyr::summarise(
      coin_results = paste(coin, collapse = "")) %>%
    dplyr::mutate(
      GT = purrr::map_chr(coin_results, function(x) {
        ifelse(paste(unique(as.vector(stringr::str_split(string = x, pattern = "", simplify = T))), collapse = "") == "H", "Ref",
               ifelse(paste(unique(as.vector(stringr::str_split(string = x, pattern = "", simplify = T))), collapse = "") == "T", "Alt",
                      "Het"))
      })
      ) %>%
    dplyr::rename(Loci = loci,
                  "Coin Results" = "coin_results",
                  "Genotype Call" = "GT") %>%
    dplyr::mutate(Locinum = as.numeric(stringr::str_split_fixed(Loci, "_", n = 2)[,2])) %>%
    dplyr::arrange(Locinum) %>%
    dplyr::select(-c("Locinum"))
}
