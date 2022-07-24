
# 00.process_DRC_DR_data.R
#
# Author: Bob Verity
# Date: 2022-07-21
#
# Purpose:
# Reads in raw data from DRC MIP analysis - taken directly from GitHib repos.
# Filters and simplifies this dataset to just a handful of mutations and saves
# simplified dataset to file.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")

# load packages
library(dplyr)
library(MIPanalyzer)

# read in data from GitHub repos associated with DRC BigBarcode data
dat <- readRDS(url("https://github.com/bobverity/antimalarial_resistance_DRC/raw/master/source_data/dr_processed2.rds"))

# subset to DRC only
dat <- filter_samples(dat, dat$samples$Country == "DRC")

# filter loci
dat <- filter_loci(dat, dat$loci$gene_name == "dhps")
dat <- filter_loci(dat, dat$loci$codon_num %in% c(437, 540, 581))

# get within-sample allele frequency of REF allele at these loci
# NB, there are no complex multi-allelic observations in this set of samples so
# we can simplify to a single frequency per locus
REF_WSAF <- matrix(NA, nrow(dat$samples), 3)
for (i in 1:nrow(REF_WSAF)) {
  x1 <- dat$counts[1,i,]
  x1[is.na(x1)] <- 0
  x2 <- dat$coverage[i,]
  x2[is.na(x2)] <- 0
  REF_WSAF[i,] <- x1 / x2
}
REF_WSAF[is.nan(REF_WSAF)] <- NA

# get trimmed version of sample dataframe
sample_trimmed <- dat$samples %>%
  select(ID, ADM1NAME)
rownames(sample_trimmed) <- NULL
names(sample_trimmed) <- c("SAMPLE_ID", "PROVINCE")

# remove samples that are NA for all 3 loci
w <- which(rowSums(is.na(REF_WSAF)) == 3)
sample_trimmed <- sample_trimmed[-w,]
REF_WSAF <- REF_WSAF[-w,]

# get trimmed version of loci dataframe
loci_trimmed <- dat$loci %>%
  select(CHROM_NUMERIC, POS, REF, ALT_list, gene_name, codon_num, codon, codon_pos)
names(loci_trimmed) <- c("CHROM", "POS", "REF", "ALT", "GENE_NAME", "CODON_NUM", "CODON", "CODON_POS")
loci_trimmed$ALT <- mapply(first, loci_trimmed$ALT)

# make final list and dataframe
dat_list <- list()
for (i in 1:nrow(sample_trimmed)) {
  dat_list[[i]] <- data.frame(SAMPLE_ID = rep(sample_trimmed$SAMPLE_ID[i], 3),
                              PROVINCE = rep(sample_trimmed$PROVINCE[i], 3)) %>%
    bind_cols(loci_trimmed) %>%
    mutate(REF_WSAF = REF_WSAF[i,])
}
dat_df <- dat_list %>%
  bind_rows()

# save to file
saveRDS(dat_df, file = "outputs/DRC_DR1.rds")
