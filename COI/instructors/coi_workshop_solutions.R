# ================================================================================
# Part 1. Loading genotype data from VCFs, exploring, and running RMCL
#
#
# ================================================================================


### Install vcfR
if (!require(vcfR)) {
    install.packages("vcfR")
}


### Set *your* working directory
setwd("/Users/jasongms/Documents/nomads/projects/AMMS2022/COI")


### Load Library
library(vcfR)


# --------------------------------------------------------------------------------
# Load VCFs and extract genotype array
#
#
# --------------------------------------------------------------------------------

# Choose files to load, and give names
vcf_files <- c(
  "data/SNPs.DRCongo.setA.random96.vcf", 
  "data/SNPs.Vietnam.setD.random96.vcf"
)
vcf_names <- c(
  "DRCongo", 
  "Vietnam"
)

# Load
vcfs <- lapply(vcf_files, read.vcfR)
names(vcfs) <- vcf_names

# Extract genotype data
gts_chr <- lapply(vcfs, extract.gt, element="GT")

# Unfortunately, GT data is as character strings; e.g. 0/0, 0/1, 1/1
# So we write a small function to convert element-wise
convert_gt_to_int <- function(gt_as_chr, divide_by=1) {
  # Handle case where the genotype is missing
  if (is.na(gt_as_chr)) {
    return(-2/divide_by)
  }
  
  # Split unphased genotype, convert to numeric, sum
  gt_as_int <- sum(as.numeric(unlist(strsplit(gt_as_chr, "/"))))
  return(gt_as_int/divide_by)
}

# Convert to ints
gts <- lapply(gts_chr, function(x) { apply(x, c(1, 2), convert_gt_to_int) })


# --------------------------------------------------------------------------------
# Plot PLAFS
#
#
# --------------------------------------------------------------------------------


calc_plafs <- function(gt_array) {
  
  # Let's make missing data NA to avoid including it in counts
  gt_missing_na = gt_array
  gt_missing_na[gt_array == -2] = NA
  
  alt_count = rowSums(gt_missing_na, na.rm=T)
  not_missing_count = 2*rowSums(!is.na(gt_missing_na)) # NB: multiply by 2, assuming samples as diploid
  
  plaf = alt_count/not_missing_count
  
  return(plaf)
}
plafs <- lapply(gts, calc_plafs)

# Plot
# - Different sets will have different distributions
hist(plafs$DRCongo, breaks=seq(0, 1, 0.05), col="Green")
hist(plafs$Vietnam, breaks=seq(0, 1, 0.05), col="Orange")


# --------------------------------------------------------------------------------
# Run THEREALMcCOIL
#
#
# --------------------------------------------------------------------------------


# Install THEREALMcCOIL
# Do this into the /COI directory
# `git clone https://github.com/EPPIcenter/THEREALMcCOIL.git`
# Also see: https://github.com/EPPIcenter/THEREALMcCOIL/tree/master/categorical_method

# Reformat data
reformat_for_rmcl <- function(gt_array) {
  
  # First make sure input format is as expected
  u_elements <- sort(unique(c(gts$Vietnam)))
  stopifnot(all(u_elements == c(-2, 0, 1, 2)))
  
  # Reformat
  gt <- gt_array / 2 # divide values by two
  gt <- t(gt)  # transpose
  
  return(gt)
}
gts_rmcl <- lapply(gts, reformat_for_rmcl)

# RUN RMCL
# Prepare output directory
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Source code, move into correct working direcory
orig_wd <- getwd()
setwd("THEREALMcCOIL/categorical_method/")
source("McCOIL_categorical.R")

# Run on each SNP set
output_paths <- c()
for (snp_set_name in names(gts_rmcl)) {
  
  # Get the array
  snp_set_gt <- gts_rmcl[[snp_set_name]]
  
  # Define output location
  # NB: you have to use a relative path, else RMCL fails
  output_path <- file.path(output_dir, paste(snp_set_name, "rmcl.txt", sep="."))
  
  # Run
  McCOIL_categorical(
    snp_set_gt,
    maxCOI=25,
    threshold_ind=20,
    threshold_site=20,
    totalrun=1000,
    burnin=100,
    M0=15,
    e1=0.05,
    e2=0.05,
    err_method=3,
    path=getwd(),
    output=file.path("..", "..", output_path) # Up from THEREALMcCOIL/categorical/
  )

  # Store path name for loading
  output_paths <- c(output_paths, output_path)
}
setwd(orig_wd)


# --------------------------------------------------------------------------------
# Explore THEREALMcCOIL outputs
#
#
# --------------------------------------------------------------------------------

# Load summary data
summary_output_paths <- sapply(output_paths, function(s) {paste(s, "summary.txt", sep="_")})
rmcl_dfs <- lapply(summary_output_paths, read.table, header=T, sep="\t")
names(rmcl_dfs) <- vcf_names

# Subset to COI results
coi_dfs <- lapply(rmcl_dfs, subset, CorP == 'C')

# Compute fraction mixed
frac_mixed <- lapply(coi_dfs, function(df) {sum(df$mean > 1) / dim(df)[1]})
barplot(unlist(frac_mixed), col=c("Green", "Orange"))

# What does the distribution look like?
hist(coi_dfs$DRCongo$mean, breaks=seq(0, 10, 1), col="Green")
hist(coi_dfs$Vietnam$mean, breaks=seq(0, 10, 1),  col="Orange")


# ... probably a lot other exploration is possible


# ================================================================================
# Part 2.Loading allelic depth data from VCF, computing WSAF, heterozygosity
# and Fws
#
# ================================================================================


# VCFs are already loaded, now we want to extract allelic depth
ads_chr <- lapply(vcfs, extract.gt, element="AD")

# Again, these are character strings: e.g. 10,45 = 10 REF reads, 45 ALT reads
# We can convert directly to WSAF
convert_ad_to_wsaf <- function(ad_as_chr, divide_by=1) {
  # Split unphased genotype, convert to numeric, sum
  ads <- as.numeric(unlist(strsplit(ad_as_chr, ",")))
  total_depth <- sum(ads)
  
  # Handle no coverage case
  if (total_depth == 0) {
    return(NA)
  }
  
  alt_depth <- ads[2] # 1-indexed
  wsaf <- alt_depth / total_depth
  
  return(wsaf)
}
wsafs <- lapply(ads_chr, function(x) { apply(x, c(1, 2), convert_ad_to_wsaf) })

# We can also compute the mean WSAF for each site
mean_wsafs <- lapply(wsafs, function(x) {rowMeans(x, na.rm=T)})

# We can compare PLAF vs. mean WSAFs
plot(x=plafs$Vietnam, y=mean_wsafs$Vietnam)
# - Linear -- not at all surprising
# - Both PLAF and mean WSAF per-site driven by # samples not REF

# Compute heterozygosities
calc_biallelic_heterozygosity <- function(p_alt) {
  return (p_alt * (1 - p_alt))
}

het_wsafs <- lapply(wsafs, calc_biallelic_heterozygosity)
het_plafs <- lapply(plafs, calc_biallelic_heterozygosity)

# Look at relationship for PLAF
plot(plafs$DRCongo, het_plafs$DRCongo)

# Finally, compute Fws
fws = lapply(vcf_names, function(vcf_name) { 
  colMeans(1 - het_wsafs[[vcf_name]] / het_plafs[[vcf_name]], na.rm=TRUE)
})
names(fws) <- vcf_names

# Visualise relationship with COI
# Add Fws column to coi_dfs
final_dfs <- lapply(vcf_names, function(vcf_name) { cbind(coi_dfs[[vcf_name]], "fws"=fws[[vcf_name]])})
names(final_dfs) <- vcf_names

# Plot
dr_final_df <- final_dfs$DRCongo
boxplot(dr_final_df$fws ~ dr_final_df$mean)
# - and we see the expected declining relationship

