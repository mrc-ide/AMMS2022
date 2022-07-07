# AMMS 2022 Tutorial: The complexities of complexity of infection (COI)
## Completed
- Created some 96 SNP sets by downsampling from MalariaGEN Pf3k data;
	- Stored /data
	- DRCongo sets have different PLAFS
	- Plus a set from Vietnam
- Some basic R code for loading VCFs, running THEREALMcCOIL, computing Fws
  - Found in `coi_workshop_solutions.R`

## TODO
- Figure out / adjust overall narrative given results
- Go from a bunch of R code into a workshop, e.g. with questions participants can move through.
  - Munging the VCF is a bit ugly, do we expect participants to have R skills to do this themselves?
  - Can we simplify?
  - How should we prepare?
	- README.md of /COI?
	- R Markdown?
	- Doc on google drive?

## Outstanding questions
1. Will participants have `git` installed? For this tutorial they will need it.
2. Is there a better interface for RMCL? E.g. an R package?
3. Should we actually use these SNP sets? I tried computing COI across different PLAF SNPs, but results looked similar; so not much learned there. Instead focusseed on DRCongo / Vietnam comparision; but maybe then we just want to use barcode SNP data. Though, we do need allelic depths to comput Fws
