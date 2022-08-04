mtdt <- readRDS("data/metadata.RDS")
mtdt <- mtdt %>%
  dplyr::group_by(deme) %>%
  dplyr::mutate(
    smpl = purrr::map_chr(1:5, function(x){paste})
  )


mtdt <- rbind(mtdt, mtdt, mtdt, mtdt, mtdt) %>%
  dplyr::arrange(deme) %>%
  dplyr::mutate(smpl = paste0(deme, 1:5)) %>%
  dplyr::select(c("deme", "smpl", dplyr::everything()))
mtdt1 <- mtdt %>%
  dplyr::select(c("deme", "longnum", "latnum")) %>%
  dplyr::filter(!duplicated(.)) %>%
  sf::st_as_sf(., coords = c("longnum", "latnum"), crs = 4326)

mtdt2 <- mtdt %>%
  dplyr::select(c("deme", "longnum", "latnum")) %>%
  dplyr::filter(!duplicated(.)) %>%
  sf::st_as_sf(., coords = c("longnum", "latnum"), crs = 4326)

gcdist <- sf::st_distance(x = mtdt1, y = mtdt2)
colnames(gcdist) <- rownames(gcdist) <- LETTERS[1:5]
gcdist %>%
  as.dist() %>%
  broom::tidy() %>%
  dplyr::rename(deme_p1 = item1,
                deme_p2 = item2) %>%
  dplyr::mutate(deme_p1 = as.character(deme_p1),
                deme_p2 = as.character(deme_p2),
                distance = distance/1e3) %>%
  saveRDS(., "data/deme_gc_dist.RDS")
