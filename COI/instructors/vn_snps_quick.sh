#!/usr/bin/env bash

bcftools concat vcfs/SNP_INDEL_Pf3D7_01_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_02_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_03_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_04_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_05_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_06_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_07_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_08_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_09_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_10_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_11_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_12_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_13_v3.combined.filtered.vcf.gz \
vcfs/SNP_INDEL_Pf3D7_14_v3.combined.filtered.vcf.gz |
vcftools --vcf - --positions chrompos.bed --keep vnsmpls.txt \
--out SNPs.Vietnam.setB.random96
