library(tidyverse)
library(data.table)
library(TwoSampleMR)
## Immune cell on GWAS Catalog Data
exp <- paste0('ebi-a-GCST9000', 1391:2121)
exp_data<- try(extract_instruments(outcomes = exp,p1 = 5e-06, clump = T, r2 = 0.001, kb = 10000))
dfOut <- fread('GCST90044134.txt.gz',data.table = F) 
out_data <- format_data(out_data,type = "outcome",snps = exp_data$SNP)
dat <- harmonise_data(exposure_dat = exp_data,outcome_dat = out_data)
dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                  dat$se.exposure,
                                  dat$samplesize.exposure)
dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
                                 dat$se.outcome,
                                 dat$samplesize.outcome)
dat <- steiger_filtering(dat)
dat$F_statistic <- round(dat$beta.exposure ^ 2 / dat$se.exposure ^ 2,2)
res <- mr(dat)
res_hete <- TwoSampleMR::mr_heterogeneity(dat)
res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
res_leaveone <- TwoSampleMR::mr_leaveoneout(dat)
res <- generate_odds_ratios(res)
res$estimate <- paste0(
  format(round(res$or,2),nsmall = 2),'(',
  format(round(res$or_lci95,2),nsmall = 2),'-',
  format(round(res$or_uci95,2),nsmall = 2),')')
rio::export(res,'All_results for GCST90044134 MR.xlsx')
res_IEU_ivw <- res %>%
  filter(., method == 'Inverse variance weighted') %>%
  filter(., pval < 0.05) %>%
  select(., exposure, nsnp, or, estimate, pval)
rio::export(res_IEU_ivw, 'IVW for GCST90044134 MR.xlsx')

## Immune cell on GWAS Catalog Data
exp <- paste0('ebi-a-GCST9000', 1391:2121)
exp_data<- try(extract_instruments(outcomes = exp,p1 = 5e-06, clump = T, r2 = 0.001, kb = 10000))
dfOut <- fread('finngen_R7_R18_DYSPHAGIA.txt.gz',data.table = F) 
out_data <- format_data(out_data,type = "outcome",snps = exp_data$SNP)
dat <- harmonise_data(exposure_dat = exp_data,outcome_dat = out_data)
dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                  dat$se.exposure,
                                  dat$samplesize.exposure)
dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
                                 dat$se.outcome,
                                 dat$samplesize.outcome)
dat <- steiger_filtering(dat)
dat$F_statistic <- round(dat$beta.exposure ^ 2 / dat$se.exposure ^ 2,2)
res <- mr(dat)
res_hete <- TwoSampleMR::mr_heterogeneity(dat)
res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
res_leaveone <- TwoSampleMR::mr_leaveoneout(dat)
res <- generate_odds_ratios(res)
res$estimate <- paste0(
  format(round(res$or,2),nsmall = 2),'(',
  format(round(res$or_lci95,2),nsmall = 2),'-',
  format(round(res$or_uci95,2),nsmall = 2),')')
rio::export(res,'All_results for finngen_R7_R18_DYSPHAGIA MR.xlsx')
res_Finngen_ivw <- res %>%
  filter(., method == 'Inverse variance weighted') %>%
  filter(., pval < 0.05) %>%
  select(., exposure, nsnp, or, estimate, pval)
rio::export(res_Finngen_ivw, 'IVW for finngen_R7_R18_DYSPHAGIA MR.xlsx')

## Take the Intersection
Intersection_IDs <- intersect(res_IEU_ivw$exposure, res_Finngen_ivw$exposure)