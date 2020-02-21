#' Get distributions of demographics
#' 

library(dplyr)

dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'
clin <- 'nationwidechildrens.org_clinical_patient_gbm.txt'
clin <- read.table(file.path(dd, clin), sep='\t', skip=1, header=T, stringsAsFactors=F) # get clinical
clin <- clin[-1,]
clin$X <- clin$bcr_patient_barcode
clin$age_cat <- factor(unlist(lapply(clin$age_at_initial_pathologic_diagnosis, function(i) ifelse(i < 50, 1, 0))))
clin$age_at_initial_pathologic_diagnosis <- as.numeric(clin$age_at_initial_pathologic_diagnosis)
clin$karnofsky_performance_score <- as.numeric(clin$karnofsky_performance_score)
clin$eastern_cancer_oncology_group <- as.numeric(clin$eastern_cancer_oncology_group)

gene <- 'gene_expression_528.txt'; gene <- read.table(file.path(dd, gene))
gene$X <- rownames(gene)
gene <- gene[, 'X', drop=F]

vasari <- 'vasari_annotations_191.csv'; vasari <- read.csv(file.path(dd, vasari), row.names=1)
vasari$X <- rownames(vasari)
vasari <- vasari[vasari$X %in% gene$X,]
nrow(vasari)

ae <- gene[!(gene$X %in% vasari$X), , drop=F]

verhaak <- 'TCGA_unified_CORE_ClaNC840.txt'; verhaak <- read.csv(file.path(dd, verhaak), sep='\t', stringsAsFactors=F)
verhaak <- as.data.frame(t(verhaak[, 3:ncol(verhaak)]))[, 1, drop=F]
attr(verhaak$subtype, 'ATT') <- NULL
colnames(verhaak) <- 'subtype'
rownames(verhaak) <- substr(rownames(verhaak), 0, 12)
rownames(verhaak) <- gsub('\\.', '\\-',rownames(verhaak))
verhaak$X <- rownames(verhaak)
nrow(verhaak)
verhaak <- verhaak[verhaak$X %in% gene$X,]
nrow(verhaak)

# combine info
clin_cols <- c('X','gender','race','ethnicity','age_at_initial_pathologic_diagnosis', 'age_cat', 'vital_status')
clin_cols <- colnames(clin)
clin_gene <- right_join(clin[, clin_cols], gene) # attach clinical info
clin_gene <- left_join(clin_gene, verhaak) # attach subtype info
clin_ae <- right_join(clin[, clin_cols], ae)
clin_ae <- left_join(clin_ae, verhaak)
clin_vasari <- right_join(clin[, clin_cols], vasari)
clin_vasari <- left_join(clin_vasari, verhaak)
clin_verhaak <- right_join(clin[, clin_cols], verhaak)

get_clin_info <- function(cohort) {
  print('gender'); print(table(cohort$gender, useNA='ifany'))
  print('race'); print(table(cohort$race, useNA='ifany'))
  print('ethnicity'); print(table(cohort$ethnicity, useNA='ifany'))
  print('age_at_initial_pathologic_diagnosis'); print( cohort %>% summarize(min=min(age_at_initial_pathologic_diagnosis,na.rm = TRUE),
                          mean=mean(age_at_initial_pathologic_diagnosis, na.rm = TRUE),
                          median=median(age_at_initial_pathologic_diagnosis, na.rm = TRUE),
                          max=max(age_at_initial_pathologic_diagnosis, na.rm = TRUE)))
  print('age_cat'); print(table(cohort$age_cat, useNA='ifany'))
  print('subtype'); print(table(cohort$subtype, useNA='ifany'))
  print('vital status'); print(table(cohort$vital_status, useNA='ifany'))
  print('initial_pathologic_diagnosis_method'); print(table(cohort$initial_pathologic_diagnosis_method, useNA='ifany'))
  # print('person_neoplasm_cancer_status'); print(table(cohort$person_neoplasm_cancer_status, useNA='ifany'))
  print('karnofsky_performance_score'); print( cohort %>% summarize(min=min(karnofsky_performance_score,na.rm = TRUE),
                                                                    mean=mean(karnofsky_performance_score, na.rm = TRUE),
                                                                    median=median(karnofsky_performance_score, na.rm = TRUE),
                                                                    max=max(karnofsky_performance_score, na.rm = TRUE)))
  print('eastern_cancer_oncology_group'); print( cohort %>% summarize(min=min(eastern_cancer_oncology_group,na.rm = TRUE),
                                                                    mean=mean(eastern_cancer_oncology_group, na.rm = TRUE),
                                                                    median=median(eastern_cancer_oncology_group, na.rm = TRUE),
                                                                    max=max(eastern_cancer_oncology_group, na.rm = TRUE)))
  # print('tissue_source_site'); print(table(cohort$tissue_source_site, useNA='ifany'))
  # check performance_status_scale_timing
}

get_clin_info(clin_gene)
get_clin_info(clin_ae)
get_clin_info(clin_vasari)
get_clin_info(clin_verhaak)

