dat <- readRDS("dr_processed2.rds")
dat.samples <- dat$samples %>% filter(Country == "DRC") %>% 
  mutate(REGION = ifelse(long<22.5,"WEST","EAST")) %>% 
  select(ID,REGION)

dat.loci <-dat$loci %>% 
  filter(gene_name == "dhps",codon_num %in% c(437, 540, 581)) %>% 
  select(POS)

dat.loci$REF_AA = c("A","K","A")
dat.loci$ALT_AA = c("G","E","G")

dr.data <- readRDS("outputs/DRC_DR1.rds")  %>% 
  left_join(dat.samples, by =c("SAMPLE_ID"="ID")) %>% 
  left_join(dat.loci , by = "POS") %>% 
  select(SAMPLE_ID,REGION,CHROM,POS,REF,ALT,GENE_NAME,CODON_NUM,CODON,CODON_POS,REF_AA,ALT_AA,REF_WSAF)

saveRDS(dr.data,"outputs/DRC_DR2.rds")

