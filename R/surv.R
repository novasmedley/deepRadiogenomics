#' Get distributions of demographics
# '

library(dplyr)
library(tidyr)
library(plyr)
library(survival)
library(foreach)
library(doParallel)
library(survminer)

library(ggplot2)
library(viridis)
library(RColorBrewer)
library(wesanderson)
library(rcartocolor)
library(scales)
library(cowplot) # arrange graphs

fdd <- '/Volumes/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_figs/clin_figs'
fdd <- file.path('~/Desktop/fig_pdfs')

#
# load data -------------------------------------------------------------------------------

dd <- '/Volumes/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'
clin <- 'nationwidechildrens.org_clinical_patient_gbm.txt'
clin <- read.table(file.path(dd, clin), sep='\t', skip=1, header=T, stringsAsFactors=F) # get clinical
clin <- clin[-1,]
clin$X <- clin$bcr_patient_barcode
clin$age_cat <- factor(unlist(lapply(clin$age_at_initial_pathologic_diagnosis, function(i) ifelse(i < 50, 1, 0))))
clin$age_at_initial_pathologic_diagnosis <- as.numeric(clin$age_at_initial_pathologic_diagnosis)
clin$karnofsky_performance_score <- as.numeric(clin$karnofsky_performance_score)
clin$eastern_cancer_oncology_group <- as.numeric(clin$eastern_cancer_oncology_group)

gene <- 'gene_expression.txt'; gene <- read.table(file.path(dd, gene))
gene$X <- rownames(gene)
gene <- gene[, 'X', drop=F]
length(gene$X[gene$X %in% clin$X]) # 502 patients have clinical info, 26 do not

vasari <- 'vasari_annotations_191.csv'; vasari <- read.csv(file.path(dd, vasari), row.names=1)
vasari$X <- rownames(vasari)
vasari <- vasari[vasari$X %in% gene$X,]
nrow(vasari)

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

#
# load label info -------------------------------------------------------------------------------

data_dd <- '/Volumes/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular/nn_retraining_data'
lab_fns <- c('f5','f6','f7','f9','f10','f14','subtype')
parallel::mcaffinity(1:7)
cl <- makeCluster(7)
registerDoParallel(cl)
labels <-  foreach (i=1:length(lab_fns), .export=c()) %dopar% {
  read.csv(file.path(data_dd, paste0(lab_fns[i],'_original.csv')), stringsAsFactors=F)[, 1:2]
}
stopCluster(cl)

labels <- Reduce(function(x,y) full_join(x, y, by='X'), labels)
colnames(labels) <- lapply(colnames(labels), function(i) strsplit(i, '\\.')[[1]][1] )

rownames(labels) <- labels$X
labels <- labels[, 2:ncol(labels)]
labels$X <- rownames(labels)


#
# TCGA-GBM info (all) -------------------------------------------------------------------------------
clin_dir <- '/Volumes/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/supplemental_materials/TCGA_clinical_data'
# drugs
drugs <- read.table(file.path(clin_dir, 'nationwidechildrens.org_clinical_drug_gbm.txt'), 
                    sep='\t', skip=1, header=T, stringsAsFactors=F) 
drugs <- drugs[-1,]
drugs$X <- drugs$bcr_patient_barcode
# drugs <- drugs[drugs$X %in% vasari$X,]
drugs <- drugs[drugs$regimen_indication %in% c('ADJUVANT'),]
length(unique(drugs$bcr_patient_barcode)) # unique patients with adjuvant drugs
table(drugs$regimen_indication)

# radiation
rad <- read.table(file.path(clin_dir, 'nationwidechildrens.org_clinical_radiation_gbm.txt'), 
                  sep='\t', skip=1, header=T, stringsAsFactors=F) 
rad <- rad[-1,]
rad$X <- rad$bcr_patient_barcode
# rad <- rad[rad$X %in% vasari$X,]
table(rad$regimen_indication)
length(unique(rad$bcr_patient_barcode)) # unique patients 

# follow-up
fu <- read.table(file.path(clin_dir, 'nationwidechildrens.org_clinical_radiation_gbm.txt'), 
                 sep='\t', skip=1, header=T, stringsAsFactors=F) 
fu <- fu[-1,]
fu$X <- fu$bcr_patient_barcode
# fu <- fu[fu$X %in% vasari$X,]
length(unique(fu$bcr_patient_barcode)) # unique patients

fu2 <- read.table(file.path(clin_dir, 'nationwidechildrens.org_clinical_follow_up_v1.0_nte_gbm.txt'), 
                  sep='\t', skip=1, header=T, stringsAsFactors=F) 
fu2 <- fu2[-1,]
fu2$X <- fu2$bcr_patient_barcode
# fu2 <- fu2[fu2$X %in% vasari$X,]
length(unique(fu2$bcr_patient_barcode)) # unique patients
table(fu2$new_neoplasm_event_type)



# -- patient events -------------------------------------------------------------------------------
# merge events data
dg <- drugs[, c('X','days_to_drug_therapy_start','days_to_drug_therapy_end', 'regimen_indication')]
rd <- rad[, c('X', 'days_to_radiation_therapy_start', 'days_to_radiation_therapy_end', 'regimen_indication')]
nte <- fu2[, c('X', 'days_to_new_tumor_event_after_initial_treatment', 
               'days_to_new_tumor_event_additional_surgery_procedure', 'new_neoplasm_event_type')]
colnames(dg) <- c('X', 'a', 'b','status'); dg$event <- 'drug'
colnames(rd) <- c('X', 'a', 'b','status'); rd$event <- 'radiation'
colnames(nte) <- c('X', 'a', 'b','status'); nte$event <- 'nte_or_surg'

dg$type <- apply(drugs[, c('drug_name'), drop=F], 1,  function(i) paste(i, collapse=':'))
dg$dose <- apply(drugs[, c('total_dose','total_dose_units')], 1, function(i) paste(i, collapse=':'))
rd$type <- rad$radiation_type
rd$dose <- rad$radiation_dosage
nte$type <- NA
nte$dose <- NA
events <- bind_rows(dg, rd, nte)

# just select last known date
events[, 2:3] <- apply(events[, 2:3], 2, function(i) as.numeric(i))
events[is.na(events)] <- -1
events$days <- pmax(events$a, events$b)
# events <- events[, c('X','days','event','status')]
events <- events %>% arrange(X, days)
events$type <- tolower(events$type)

#
# -- km time to death, get info -------------------------------------------------------------------------------

km <- clin[, c('X','vital_status', 'days_to_last_followup','days_to_death')]
km <- km %>% 
  mutate(outcome = ifelse(vital_status == "Alive", 0, 1))

km <- split(km, km$vital_status)
km[['Dead']]$days <- as.numeric(km[['Dead']]$days_to_death)
km[['Dead']] <- km[['Dead']][, c('X','days','outcome')]

# update survival days by checking follow up days in other metadata files
temp <- km[['Alive']][c('X','vital_status','days_to_last_followup')]
colnames(temp) <- c('X','status','days')
temp$event <- 'initial'
temp$days <- as.numeric(temp$days)
temp <- bind_rows(temp, events[events$X %in% temp$X, ]) # look at their other events
temp <- temp %>% group_by(X) %>%
  arrange(desc(days)) %>%
  slice(1L) # get the last followup day for each patient

# merge alive info and dead info
temp$days <- as.numeric(temp$days)
temp$outcome <- 0 # all were alive patients
km <- bind_rows(temp[, colnames(temp) %in% colnames(km[['Dead']])], km[['Dead']])
km$yrs <- km$days/365 # convert to years

km <- km[km$X %in% gene$X, ]
#
# -- km time to death, plot -------------------------------------------------------------------------------
dev.new() 

# TCGA-GBM survival
fit_all_os <- survfit(Surv(km$yrs, km$outcome) ~ 1); fit_all_os
plot(fit_all_os, xlab = "Years", ylab = "Overall survival probability")  
ggsurvplot(fit=fit_all_os, xlab = "Years", ylab = "Overall survival probability",
           data = km, break.x.by = 1, censor = FALSE, risk.table = TRUE, risk.table.y.text = FALSE)

# TCGA-GBM, VASARI survival
vaskm <- bind_cols(vasari[vasari$X %in% km$X, ], km[km$X %in% vasari$X, ])
fit_vasari_os <- survfit(Surv(vaskm$yrs, vaskm$outcome) ~ 1); fit_vasari_os
dev.new(); plot(fit_vasari_os, xlab = "Years", ylab = "Overall survival probability")  
ggsurvplot(fit=fit_vasari_os, xlab = "Years", ylab = "Overall survival probability",
           data = vaskm, break.x.by = 1, censor = FALSE, risk.table = TRUE, risk.table.y.text = FALSE)

# per trait splits
get_km_s <- function(df, y_colname, title) {
  form <- as.formula(paste('Surv(yrs, outcome) ~', y_colname))
  f <- do.call(survfit, list(formula=form, data=df))
  print(ggsurvplot(fit=f, xlab = "Years", ylab = "Overall survival probability", legend.title = title,
                   data = df, break.x.by = 1, censor = FALSE, risk.table = TRUE, risk.table.y.text = FALSE))
  dif <- do.call(survdiff, list(formula=Surv(df$yrs, df$outcome) ~ df[, y_colname], data=df))
  # print(dif)
  
  p.val <- 1 - pchisq(dif$chisq, length(dif$n) - 1)
  print(p.val)
  
  f$cohort <- y_colname
  f$p.val <- p.val
  return(f)
}
dev.new();
cn <- 'f5'; fit_f5_os <- get_km_s(df=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), y_colname=cn, cn)
cn <- 'f6'; fit_f6_os <-get_km_s(df=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), y_colname=cn, cn)
cn <- 'f7'; fit_f7_os <-get_km_s(df=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), y_colname=cn, cn)
cn <- 'f14'; fit_f14_os <-get_km_s(df=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), y_colname=cn, cn)
cn <- 'f9'; fit_f9_os <-get_km_s(df=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), y_colname=cn, cn)
cn <- 'f10'; fit_f10_os <-get_km_s(df=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), y_colname=cn, cn)

#
# -- km time to new tumor event, get info -------------------------------------------------------------------------------

# get patients with progression
prog <- events[events$event=='nte_or_surg', c('X','status','days')]
prog[prog==-1] <- NA
prog <- na.omit(prog)
table(prog$status)
unknown <- prog[prog$status=='[Unknown]',] 

#' 6 unknowns in vasari dataset (after checking events, prog, and km data):
#' 0176 has other prog events, has progression --> remove event, 1
#' 0184 has unknown at 1276, but radiation at 1300 with progressions status, has progression --> set unknown as progression, 1
#' 0188 has unknown, no other event info expectp alive in km data --> remove event, 0
#' 0241 has recurrence on same day as unknown --> remove event, 1
#' 0644 ""
#' 1802 ""
#' 
#' 21 unknowns not in vasari:
#' 0744 no other data --> remove event, 0
#' 0876 has unknown at 424, but radiation at 446 with progression status --> set unknown as progression, 1
#' 0939 has progression --> remove event, 1
#' 1084 ""
#' 1801 ""
#' 2563 has unknown at 554, but radiation at 576 with progression status --> set unknown as progression, 1
#' 5956 has progression --> remove event, 1
#' 5208 ""
#' 5211 ""
#' 5215 ''
#' 5219 ''
#' 5220 ''
#' 4213 ''
#' 5222 ''
#' 5651 ''
#' 6575 ''
#' 6577 ''
#' 6578 ''
#' 6584 ''

# set unknowns with radiation inform stating progression
prog$status[prog$X %in% c('TCGA-06-0184', 'TCGA-06-0876', 'TCGA-06-2563')] <- 'Progression'
nrow(unknown)
unknown <- prog[prog$status=='[Unknown]',] 
nrow(unknown)

prog <- prog[prog$status!='[Unknown]',] # remove unknown
prog <- prog %>% arrange(X, days)
prog <- prog %>% group_by(X) %>%
  arrange(days) %>%
  slice(1L) # get day of first biopsy for each patient
table(prog$status)
length(unique(prog$X)); nrow(prog)

# patients without progression??
no_prog <- km[!(km$X %in% prog$X), ]
no_prog <- no_prog[, c('X','days')]
no_prog$outcome <- 0
length(unique(no_prog$X))

# combine prog with no prog patients
prog$outcome <- 1
prog <- bind_rows(prog, no_prog)
prog <- prog[, c('X', 'days','outcome')]
prog$yrs <- prog$days/365
prog <- prog[prog$X %in% gene$X, ]

#
# -- km time to new tumor event, plot -------------------------------------------------------------------------------
dev.new(); 

# TCGA-GBM progression
fit_all_prog <- survfit(Surv(prog$yrs, prog$outcome) ~ 1); fit_all_prog
fit_all_prog$cohort <- 'TCGA-GBM'
fit_all_prog$p.val <- 'NA'
plot(fit_all_prog, xlab = "Years", ylab = "Overall progression free survival probability")  
# ggsurvplot(fit=fit, xlab = "Years", ylab = "Overall progression  free survival probability",
#            data = fit_all_prog, break.x.by = 1, censor = FALSE, risk.table = TRUE, risk.table.y.text = FALSE)

# TCGA-GBM, VASARI progression
vasprog <- bind_cols(vasari[vasari$X %in% prog$X, ], prog[prog$X %in% vasari$X, ])
fit_vasari_prog <- survfit(Surv(vasprog$yrs, vasprog$outcome) ~ 1); fit_vasari_prog
dev.new(); plot(fit_vasari_prog, xlab = "Years", ylab = "Overall progression  free survival probability")  
ggsurvplot(fit=fit_vasari_prog, xlab = "Years", ylab = "Overall progression free survival probability",
           data = vasprog, break.x.by = 1, censor = FALSE, risk.table = TRUE, risk.table.y.text = FALSE)

# per trait splits progression
dev.new();

get_km_p <- function(df, y_colname, title) {
  form <- as.formula(paste('Surv(yrs, outcome) ~', y_colname))
  f <- do.call(survfit, list(formula=form, data=df))
  print(ggsurvplot(fit=f, xlab = "Years", ylab = "Overall progression free survival probability", legend.title = title,
                   data = df, break.x.by = 1, censor = FALSE, risk.table = TRUE, risk.table.y.text = FALSE))
  dif <- do.call(survdiff, list(formula=Surv(df$yrs, df$outcome) ~ df[, y_colname], data=df))
  # print(dif)
  
  p.val <- 1 - pchisq(dif$chisq, length(dif$n) - 1)
  print(p.val)
  
  f$cohort <- y_colname
  f$p.val <- p.val
  return(f)
}

cn <- 'f5'; fit_f5_prog <- get_km_p(df=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), y_colname=cn, cn)
cn <- 'f6'; fit_f6_prog <- get_km_p(df=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), y_colname=cn, cn)
cn <- 'f7'; fit_f7_prog <- get_km_p(df=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), y_colname=cn, cn)
cn <- 'f14'; fit_f14_prog <- get_km_p(df=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), y_colname=cn, cn)
cn <- 'f9'; fit_f9_prog <- get_km_p(df=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), y_colname=cn, cn)
cn <- 'f10'; fit_f10_prog <- get_km_p(df=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), y_colname=cn, cn)




#
# -- cox prep -------------------------------------------------------------------------------

# get non-imaging variables
diag_median_age <- mean(clin$age_at_initial_pathologic_diagnosis, na.rm = TRUE)
these <-  c('X','gender','race','ethnicity','age_at_initial_pathologic_diagnosis', 'initial_pathologic_diagnosis_method')
clinv <- clin[clin$X %in% vasari$X, these]
clinv$age_at_initial_pathologic_diagnosis <- ifelse(clinv$age_at_initial_pathologic_diagnosis < diag_median_age, 'BELOW', 'ABOVE')
clinv <- left_join(clinv, verhaak)

# clean up
clinv$subtype <- as.character(clinv$subtype)
# clinv[is.na(clinv)] <- '[Not Available]'
table(clinv$subtype)

clinv$race <- mapvalues(clinv$race, from=c('ASIAN','BLACK OR AFRICAN AMERICAN', '[Not Available]'),
                        to=c('NON-WHITE', 'NON-WHITE', NA))
table(clinv$race)
clinv$ethnicity <- mapvalues(clinv$ethnicity, from=c('HISPANIC OR LATINO','NOT HISPANIC OR LATINO', '[Not Available]'),
                             to=c('HISPANIC', 'NON-HISPANIC', NA))
table(clinv$ethnicity)

clinv$initial_pathologic_diagnosis_method <- mapvalues(clinv$initial_pathologic_diagnosis_method,
                                                       from=c('Excisional Biopsy', 'Fine needle aspiration biopsy',
                                                              'Other method, specify:', '[Not Available]',
                                                              'Tumor resection'),
                                                       to=c('OTHER', 'OTHER', 'OTHER', '[Not Available]', 'RESECTION'))
table(clinv$initial_pathologic_diagnosis_method)



# add labels info
clinv <- left_join(clinv, labels[, !(colnames(labels) %in% c('subtype'))])
clinv <- left_join(clinv, vaskm[, c('X', 'yrs','outcome')])
clinv <- left_join(clinv, vasprog[, c('X', 'yrs','outcome')])
rownames(clinv) <- clinv$X
clinv <- clinv[, !(colnames(clinv) %in% 'X')] # remove ids
colnames(clinv)[1:6] <- c('gender','race','ethnicity','diag_age','diag_meth','subtype') # rename vars
clinv[1:12] <- lapply(colnames(clinv)[1:12], function(i) as.factor(clinv[, i]))

#
# -- cox model -------------------------------------------------------------------------------

# convert categorical to binary
prepcox <- function(type, remove, form, df){
  # df needs rownames
  outinfo <- c('X','yrs', 'outcome')
  if (type == 'surv') {
    outcome <- vaskm[, outinfo]
  } else {
    outcome <- vasprog[, outinfo]
  }
  z <- model.matrix(form, df[, !(colnames(df) %in% remove)])
  colnames(z)
  
  # add back labels
  z <- as.data.frame(z)
  z$X <- rownames(z)
  z <- left_join(z, outcome)
  z <- z[, !(colnames(z) %in% 'X')] # remove ids
  return(z)
}

remove <- c('X')

# use all img traits
form <- as.formula(paste0('outcome ~ gender + race + diag_age + f5 + f6 + f7 + f9 + f10 + f14')); print(form)

z <- prepcox(type='surv', remove, form, clinv) # survival
os <- coxph(Surv(yrs, outcome) ~ . , data=z)
summary(os)
os_step <- step(coxph(Surv(yrs, outcome) ~ . , data=z))
summary(os_step) # Concordance= 0.517 

z <- prepcox(type='prog', remove, form, clinv) # progression
pr <- coxph(Surv(yrs, outcome) ~ . , data=z)
summary(pr)
pr_step <- step(coxph(Surv(yrs, outcome) ~ . , data=z))
summary(pr_step) # Concordance= 0.522

# univariate
get_cox_uni <- function(label, type) {
  form <- as.formula(paste0('outcome ~ gender + race + diag_age +', label)); print(form)
  z <- prepcox(type, remove, form, clinv) # survival
  os <- coxph(Surv(yrs, outcome) ~ . , data=z)
  summary(os)
}
get_cox_uni('f5', type='surv')
get_cox_uni('f6', type='surv')
get_cox_uni('f7', type='surv')
get_cox_uni('f9', type='surv')
get_cox_uni('f10', type='surv')
get_cox_uni('f14', type='surv')

get_cox_uni('f5', type='prog')
get_cox_uni('f6', type='prog')
get_cox_uni('f7', type='prog')
get_cox_uni('f9', type='prog')
get_cox_uni('f10', type='prog')
get_cox_uni('f14', type='prog')

#
# -- radiogenomics -------------------------------------------------------------------------------
library(foreach)
library(doParallel)

dd <- '/Volumes/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular'
fdd <- '/Volumes/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_figs/gsea'
fdd <- file.path('~/Desktop/fig_pdfs')
gs_dd <- '/Volumes/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'

get_counts <- function(gs_name, count_thres=10) {
  vasari_labs <- c('f5','f6','f7','f9','f10','f14')
  counts <- lapply(vasari_labs, function(lab) {
    f <- file.path(dd,'gsea', paste('saliency', lab, gs_name, '.rds', sep='_')); print(f)
    s <- readRDS(f)
    
    c <- s$sal
    c$enriched <- ifelse(c$padj < 0.05, 1, 0) 
    
    if (length(unique(c$enriched)) > 1) {
      # get paths with enrichements greater than threshold
      path_counts <- c %>% dplyr::group_by(pathway) %>% 
        dplyr::summarize(count=sum(enriched)) %>%
        dplyr::select(pathway, count)
      c <- c[c$pathway %in% path_counts$pathway[path_counts$count >= count_thres], ]
      c <- c[, c('pathway','name','enriched')]
      colnames(c) <- c('pathway','X','enriched')
      c$X <- gsub('\\.', '-', c$X)
      c <- c %>% spread(pathway, enriched)
    } else {
      print('no significant associations')
      c <- NULL
    }
    c
  })
  names(counts) <- vasari_labs
  return(counts)
}



get_cox_uni_rg <- function(label, type, rg, remove=c('X'), covars='gender + race + diag_age +', step=F) {
  # cox prep for radiogenomic vars
  clinv$X <- rownames(clinv)
  
  pathnames <- colnames(rg)[colnames(rg) != 'X']
  if (length(pathnames)<1){
    pathnames <- ''
  } else {
    pathnames <- paste('+', paste(pathnames, collapse=' + '))
  }
  rg <- left_join(rg, clinv)
  rownames(rg) <- rg$X
  form <- as.formula(paste0('outcome ~ ', covars, label, pathnames )); print(form)
  z <- prepcox(type, remove, form, rg) # survival
  
  if (step) {
    nothing <- coxph(Surv(yrs, outcome) ~ genderMALE + raceWHITE + diag_ageBELOW, data=z)
    full <-  coxph(Surv(yrs, outcome) ~ . , data=z)
    os <- step(nothing, scope=c(lower=as.formula(nothing), upper=as.formula(full)))
  } else {
    os <- coxph(Surv(yrs, outcome) ~ . , data=z)
  }
  print(summary(os))
}


gs_name <- 'hallmark' 
count_thres <- 5
hallmark <- get_counts(gs_name, count_thres); x <- hallmark

# cox rg
label <- 'f6'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) # check
label <- 'f7'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) # check
label <- 'f10'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) 
label <- 'f14'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) # sig vars

# survival 
cn <- 'f10'; df <- x[[cn]]; colnames(df)
gs<- 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'; z <- get_km_s(df=left_join(df[df$X %in% vaskm$X, c('X','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')], vaskm),  
                                                                 y_colname=gs, title=paste(cn, ' - ', gs))

cn <- 'f14'; df <- x[[cn]]; colnames(df)
gs<- 'HALLMARK_E2F_TARGETS'; z <- get_km_s(df=left_join(df[df$X %in% vaskm$X, c('X','HALLMARK_E2F_TARGETS')], vaskm), y_colname=gs, title=paste(cn, ' - ', gs))
gs<- 'HALLMARK_G2M_CHECKPOINT'; z <- get_km_s(df=left_join(df[df$X %in% vaskm$X, c('X','HALLMARK_G2M_CHECKPOINT')], vaskm), y_colname=gs, title=paste(cn, ' - ', gs))

# prog
cn <- 'f6'; df <- x[[cn]]; colnames(df)
gs<- 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'
z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')], vasprog), 
              y_colname=gs, title=paste(cn, ' - ', gs))
gs<- 'HALLMARK_GLYCOLYSIS'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','HALLMARK_GLYCOLYSIS')], vasprog), 
                                          y_colname=gs, title=paste(cn, ' - ', gs))
gs<- 'HALLMARK_MYOGENESIS'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','HALLMARK_MYOGENESIS')], vasprog), 
                                          y_colname=gs, title=paste(cn, ' - ', gs))

cn <- 'f7'; df <- x[[cn]]; colnames(df)
gs <- 'HALLMARK_INTERFERON_ALPHA_RESPONSE'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','HALLMARK_INTERFERON_ALPHA_RESPONSE')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs))
gs <- 'HALLMARK_P53_PATHWAY'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','HALLMARK_P53_PATHWAY')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs))

cn <- 'f10'; df <- x[[cn]]; colnames(df)
gs <- 'HALLMARK_TNFA_SIGNALING_VIA_NFKB'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','HALLMARK_TNFA_SIGNALING_VIA_NFKB')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs))





gs_name <- 'puchalski' 
count_thres <- 5
puchalski <- get_counts(gs_name, count_thres); x <- puchalski

# cox rg
label <- 'f6'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])
label <- 'f7'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) # check
label <- 'f9'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])
label <- 'f10'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])
label <- 'f14'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) # sig

# prog
cn <- 'f14'; df <- x[[cn]]; colnames(df)
gs<- 'DarmanisBarres_FetalNeuronsreplicating'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','DarmanisBarres_FetalNeuronsreplicating')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs))

gs_name <- 'verhaak' 
count_thres <- 5
verhaak_counts <- get_counts(gs_name, count_thres); x <- verhaak_counts

label <- 'f6'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])
label <- 'f7'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])
label <- 'f10'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])
label <- 'f14'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])

# prog
cn <- 'f6'; df <- x[[cn]]; colnames(df)
gs<- 'CL'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','CL')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs))
gs<- 'NL'; fit_f6_prog_NL <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','NL')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs)) # sig

cn <- 'f7'; df <- x[[cn]]; colnames(df)
gs<- 'CL'; z <- get_km_p(df=left_join(df[df$X %in% vasprog$X, c('X','CL')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs))


gs_name <- 'chromosome' 
count_thres <- 5
chromosome <- get_counts(gs_name, count_thres); x <- chromosome

label <- 'f6'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) # sig
label <- 'f7'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]])
label <- 'f10'; get_cox_uni_rg(label=label, type='surv', rg=x[[label]]) # sig

label <- 'f6'; get_cox_uni_rg(label=label, type='prog', rg=x[[label]])
label <- 'f7'; get_cox_uni_rg(label=label, type='prog', rg=x[[label]])
label <- 'f10'; get_cox_uni_rg(label=label, type='prog', rg=x[[label]])

# survival
cn <- 'f6'; df <- x[[cn]]; colnames(df)
gs<- 'chr3p21'; fit_f6_surv_chr3p21 <- get_km_s(df=left_join(df[df$X %in% vaskm$X, c('X','chr3p21')], vaskm), y_colname=gs, title=paste(cn, ' - ', gs))  #sig

cn <- 'f10'; df <- x[[cn]]; colnames(df)
gs <- 'chr19p13'; fit_f10_chr19p13 <- get_km_s(df=left_join(df[df$X %in% vaskm$X, c('X','chr19p13')], vaskm), y_colname=gs, title=paste(cn, ' - ', gs))  #sig
gs <- 'chr1p35'; fit_f10_chr1p35 <- get_km_s(df=left_join(df[df$X %in% vaskm$X, c('X','chr1p35')], vaskm), y_colname=gs, title=paste(cn, ' - ', gs))  


gs_name <- 'reactome' 
count_thres <- 5
reactome <- get_counts(gs_name, count_thres); x <- reactome

cn <- 'f5'; df <- x[[cn]]; colnames(df)
gs <- 'REACTOME_TRAFFICKING_OF_AMPA_RECEPTORS'; 
z <- get_km_s(df=left_join(df[df$X %in% vaskm$X, c('X','REACTOME_TRAFFICKING_OF_AMPA_RECEPTORS')], vaskm), y_colname=gs, title=paste(cn, ' - ', gs))  #sig
z <- get_km_s(df=left_join(df[df$X %in% vasprog$X, c('X','REACTOME_TRAFFICKING_OF_AMPA_RECEPTORS')], vasprog), y_colname=gs, title=paste(cn, ' - ', gs))  #sig

#
# COMBINE RADIOGENOMIC FEATURES ---------
# feature selection
merge_data <- function(gene_set) {
  gs <- lapply(names(gene_set), function(i) {
    colnames(gene_set[[i]]) <- paste0(i, '.', colnames(gene_set[[i]]))
    colnames(gene_set[[i]])[1] <- 'X'
    gene_set[[i]]
  })
  gs <- gs %>%
    Reduce(function(x,y) left_join(x,y, by='X'), .)
  gs
}

chr <- merge_data(chromosome[c('f6','f7', 'f10','f14')])
puch <- merge_data(puchalski[c('f6','f7', 'f9','f10','f14')])
hall <- merge_data(hallmark[c('f6','f7', 'f10','f14')])
react <- merge_data(reactome[c('f5')]) # not signifcant

df <- list(chr, puch, hall, react) %>%
  Reduce(function(x,y) inner_join(x,y, by='X'), .)
df <- na.omit(df)

rg_counts <- gather(df[, 2:ncol(df)]) %>% group_by(value) %>% dplyr::count(key)
rg_counts$key <- tolower(rg_counts$key )
rg_counts$key <- gsub('.chr','.chromosome chr', rg_counts$key)
rg_counts$key <- gsub('.darmanisbarres_','.types Darminis_', rg_counts$key)
rg_counts$key <- gsub('.barres_','.types Zhang_', rg_counts$key)
rg_counts$key <- gsub('.patel_','.types Patel_', rg_counts$key)
rg_counts$key <- gsub('.hallmark_','.hallmark ', rg_counts$key)
rg_counts$key <- gsub('.reactome_','.reactome ', rg_counts$key)
rg_counts$key <- gsub('f5','enhancing', rg_counts$key)
rg_counts$key <- gsub('f6','nCET', rg_counts$key)
rg_counts$key <- gsub('f7','necrosis', rg_counts$key)
rg_counts$key <- gsub('f9','focal', rg_counts$key)
rg_counts$key <- gsub('f10','infiltrative', rg_counts$key)
rg_counts$key <- gsub('f14','edema', rg_counts$key)

a <- strsplit(rg_counts$key, '\\.')
rg_counts$img <- unlist(lapply(a, '[[', 1))
a <- lapply(a, function(i) i[2])
a <- strsplit(unlist(a), ' ')
rg_counts$gene_set <- unlist(lapply(a, function(i) i[2]))
rg_counts$gene_set <- gsub('_',' ', rg_counts$gene_set)
rg_counts$collection <- unlist(lapply(a, function(i) i[1]))
rg_counts$collection <- gsub('type','cell types or phenotypes', rg_counts$collection)
rg_counts$enriched <- as.factor(mapvalues(as.character(rg_counts$value), from=c('0', '1'), to=c('no', 'yes')))
rg_counts <- rg_counts[rg_counts$value==1,]
# theme
t <- 10; t2 <- t + 0
th <- theme(legend.text=element_text(size=t2, color="grey30"),
            legend.title=element_text(size=t),
            legend.position='right',
            legend.key.size=unit(.8,'line'),
            legend.margin=margin(l = -0.15, unit='cm'),
            legend.background = element_rect(fill = "transparent", colour = "transparent"),
            axis.text.y=element_text(size=t2, color="#666666"),
            axis.text.x=element_text(size=t2,color='#666666', angle=90, hjust=1),
            axis.title.y = element_text(margin=margin(0,5,0,0),size=t, angle=0, vjust=0.5),
            # axis.title.x = element_text(margin=margin(5,0,0,0),size=t),
            axis.title.x = element_blank(),
            strip.text.x = element_text(color='white',size=t2,margin = margin(.05,0,.05,0, "cm")),
            strip.text.y = element_text(color='white',size=t2,margin = margin(.05,0,.05,0, "cm")),
            strip.background =  element_rect(fill='black'),
            panel.border = element_rect(colour = "black", fill=NA, size=.3),
            # panel.grid = element_line(color='grey60'),
            panel.grid.minor = element_line(colour="grey90", size=.25),
            panel.grid.major.x = element_blank()) 

ggplot(rg_counts, aes(x=gene_set, y=n, fill=enriched)) +
  geom_bar(stat='identity', position=position_dodge2(width=0.9, preserve='single', padding=0), color = "black", lwd=0.3) +
  theme_bw() + th + 
  # scale_y_continuous(limits=c(0,60), breaks=seq(0,60, 10)) +  
  scale_fill_manual(values = wes_palette("Royal1")[2]) +
  guides(fill=F) +
  facet_grid(img~collection, scale='free_x')

ff <- file.path(fdd, 'radiogenomic_traits.pdf'); print(ff)
dev.print(pdf, ff, height=9, width=8)

# survival
set.seed(425); get_cox_uni_rg(label='f5+f6+f7+f9+f10+f14', type='surv', rg=df, step=T) # radiogenomic
set.seed(425); get_cox_uni_rg(label='f5+f6+f7+f9+f10+f14', type='prog', rg=df, step=T) # radiogenomic

# univariate with same group
# survival
clinv$X <- rownames(clinv)
dfs <- left_join(df, clinv[, c('X','gender','race','diag_age','f10', 'f14', 'f6','f7', 'f5','f9')]) 
dfs <- na.omit(dfs)
get_cox_uni_rg(label='gender', type='surv', rg=dfs[, c('X'), drop=F], covars='') 
get_cox_uni_rg(label='race', type='surv', rg=dfs[, c('X'), drop=F], covars='')
get_cox_uni_rg(label='diag_age', type='surv', rg=dfs[, c('X'), drop=F], covars='')
get_cox_uni_rg(label='f10', type='surv', rg=dfs[, c('X'), drop=F], covars='') 

get_cox_uni_rg(label='f10.chr1p35',    type='surv', rg=dfs[, c('X', 'f10.chr1p35'), drop=F], covars='') 
get_cox_uni_rg(label='f14.Barres_endothelial', type='surv', rg=dfs[, c('X', 'f14.Barres_endothelial'), drop=F], covars='')
get_cox_uni_rg(label='f7.Barres_GBMcoreastrocytes', type='surv', rg=dfs[, c('X', 'f7.Barres_GBMcoreastrocytes'), drop=F], covars='')
get_cox_uni_rg(label='f7.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', type='surv', rg=dfs[, c('X', 'f7.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'), drop=F], covars='')
get_cox_uni_rg(label='f6.HALLMARK_MYOGENESIS', type='surv', rg=dfs[, c('X', 'f6.HALLMARK_MYOGENESIS'), drop=F], covars='')
get_cox_uni_rg(label='f7.HALLMARK_MYC_TARGETS_V2', type='surv', rg=dfs[, c('X', 'f7.HALLMARK_MYC_TARGETS_V2'), drop=F], covars='')
get_cox_uni_rg(label='f10.HALLMARK_MTORC1_SIGNALING', type='surv', rg=dfs[, c('X', 'f10.HALLMARK_MTORC1_SIGNALING'), drop=F], covars='')

# progression
get_cox_uni_rg(label='gender', type='prog', rg=dfs[, c('X'), drop=F], covars='')  # univariate with same group
get_cox_uni_rg(label='race', type='prog', rg=dfs[, c('X'), drop=F], covars='')
get_cox_uni_rg(label='diag_age', type='prog', rg=dfs[, c('X'), drop=F], covars='')
get_cox_uni_rg(label='f10', type='prog', rg=dfs[, c('X'), drop=F], covars='') 

get_cox_uni_rg(label='f10.chr1p35',    type='prog', rg=dfs[, c('X', 'f10.chr1p35'), drop=F], covars='') 
get_cox_uni_rg(label='f10.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',    type='prog', rg=dfs[, c('X', 'f10.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'), drop=F], covars='') 
get_cox_uni_rg(label='f7.HALLMARK_MYC_TARGETS_V2',    type='prog', rg=dfs[, c('X', 'f7.HALLMARK_MYC_TARGETS_V2'), drop=F], covars='') 
get_cox_uni_rg(label='f14.DarmanisBarres_FetalNeuronsreplicating',    type='prog', rg=dfs[, c('X', 'f14.DarmanisBarres_FetalNeuronsreplicating'), drop=F], covars='') 
get_cox_uni_rg(label='f10.HALLMARK_TGF_BETA_SIGNALING',    type='prog', rg=dfs[, c('X', 'f10.HALLMARK_TGF_BETA_SIGNALING'), drop=F], covars='') 
get_cox_uni_rg(label='f14.chr18p11',    type='prog', rg=dfs[, c('X', 'f14.chr18p11'), drop=F], covars='') 
get_cox_uni_rg(label='f14.HALLMARK_G2M_CHECKPOINT',    type='prog', rg=dfs[, c('X', 'f14.HALLMARK_G2M_CHECKPOINT'), drop=F], covars='') 
get_cox_uni_rg(label='f6.chr22q13',    type='prog', rg=dfs[, c('X', 'f6.chr22q13'), drop=F], covars='') 
get_cox_uni_rg(label='f10.chr6q27',    type='prog', rg=dfs[, c('X', 'f10.chr6q27'), drop=F], covars='') 
get_cox_uni_rg(label='f7.HALLMARK_P53_PATHWAY',    type='prog', rg=dfs[, c('X', 'f7.HALLMARK_P53_PATHWAY'), drop=F], covars='') 

# sanity check of the selected features:
temp <- dfs[, c('gender','race','diag_age','f10',
                'f10.chr1p35','f10.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                'f7.HALLMARK_MYC_TARGETS_V2','f14.DarmanisBarres_FetalNeuronsreplicating',
                'f10.HALLMARK_TGF_BETA_SIGNALING','f14.chr18p11','f14.HALLMARK_G2M_CHECKPOINT','f6.chr22q13','f10.chr6q27','f7.HALLMARK_P53_PATHWAY','X')]
get_cox_uni_rg(label='f10 +f10.chr1p35 +f10.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION +f7.HALLMARK_MYC_TARGETS_V2 +f14.DarmanisBarres_FetalNeuronsreplicating +f10.HALLMARK_TGF_BETA_SIGNALING +f14.chr18p11 +f14.HALLMARK_G2M_CHECKPOINT +f6.chr22q13 +f10.chr6q27 +f7.HALLMARK_P53_PATHWAY', type='prog', rg=temp, step=F)

# check trait freq in cox models
table(dfs$gender)
table(dfs$race)

table(dfs$f10)
table(dfs$f10.chr1p35)
table(dfs$f10.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)
table(dfs$f7.HALLMARK_MYC_TARGETS_V2)
table(dfs$f14.DarmanisBarres_FetalNeuronsreplicating)
table(dfs$f10.HALLMARK_TGF_BETA_SIGNALING)
table(dfs$f14.chr18p11)
table(dfs$f14.HALLMARK_G2M_CHECKPOINT)
table(dfs$diag_age)
table(dfs$diag_age)
table(dfs$diag_age)

#
# COMBINE PLOTS ----------------------------

#
# survival use ggsurvplot----------
add_theme <- theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
                   legend.text=element_text(size=11, color="grey30"),
                   legend.title=element_text(size=11),
                   legend.margin=margin(l = -0.15, unit='cm'),
                   legend.background = element_rect(fill = "transparent", colour = "transparent"),
                   axis.title.y = element_text(margin=margin(0,5,0,0),size=10, angle=0, vjust=0.5))

get_plot_km <- function(fit, data, title, leg.labs, pal, ylab='OS\nprobability', xlim=c(0,6)) {
  p <- ggsurvplot(fit=fit, data=data, # SURVIVAL DATA
                  xlab="years", ylab=ylab, legend.labs=leg.labs, palette=pal, xlim=xlim,
                  size=.75, break.x.by=1, censor=T, legend.title=title, legend=c(.6,.9), 
                  pval=T, pval.size=4, pval.coord=c(4.05, .55), conf.int=F,
                  risk.table=T, risk.table.y.text=F, risk.table.height=0.25, 
                  ggtheme=theme_bw(), tables.theme=theme_bw(), fontsize=3.5)
  p$table <- ggpubr::ggpar(p$table, font.title=list(size=10) ) + add_theme 
  p$plot <- p$plot + add_theme + guides(color=guide_legend(override.aes=list(size=2.5)))
  p
}
dev.new();


#
# COMBINE imaging kms ---------

# surv
cn <- 'f5'; names(fit_f5_os$strat) # must match legend.labs
plot_f5_os <- get_plot_km(fit=fit_f5_os, data=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), title='               ', 
                          leg.labs=c('enhancing < 1/3', 'enhancing >= 1/3'),  pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_f5_os)

cn <- 'f6'; names(fit_f6_os$strat) 
plot_f6_os <- get_plot_km(fit=fit_f6_os, data=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), title='               ', 
                          leg.labs=c('nCET < 1/3', 'nCET >= 1/3'),  pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_f6_os)

cn <- 'f7'; names(fit_f7_os$strat) 
plot_f7_os <- get_plot_km(fit=fit_f7_os, data=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), title='               ', 
                         leg.labs=c('necrosis < 1/3', 'necrosis >= 1/3'),  pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_f7_os)

cn <- 'f9'; names(fit_f9_os$strat) 
plot_f9_os <- get_plot_km(fit=fit_f9_os, data=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), title='               ', 
                         leg.labs=c('focal', 'non-focal'),  pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_f9_os)

cn <- 'f10'; names(fit_f10_os$strata) 
plot_f10_os <- get_plot_km(fit=fit_f10_os, data=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), title='               ', 
                           leg.labs=c('expansive', 'infiltrative'), pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_f10_os)

cn <- 'f14'; names(fit_f14_os$strat) 
plot_f14_os <- get_plot_km(fit=fit_f14_os, data=bind_cols(labels[labels$X %in% vaskm$X, c('X',cn)], vaskm), title='               ', 
                         leg.labs=c('edema < 1/3', 'edema >= 1/3'),  pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_f14_os)

a <- plot_grid(plot_f5_os$plot, plot_f5_os$table, nrow=2, rel_heights=c(1, .4))
b <- plot_grid(plot_f6_os$plot, plot_f6_os$table, nrow=2, rel_heights=c(1, .4))
c <- plot_grid(plot_f7_os$plot, plot_f7_os$table, nrow=2, rel_heights=c(1, .4))
d <- plot_grid(plot_f9_os$plot, plot_f9_os$table, nrow=2, rel_heights=c(1, .4))
e <- plot_grid(plot_f10_os$plot, plot_f10_os$table, nrow=2, rel_heights=c(1, .4))
f <- plot_grid(plot_f14_os$plot, plot_f14_os$table, nrow=2, rel_heights=c(1, .4))

img_os <- plot_grid(a, b, c, d, e, f, ncol=1)


# prog
cn <- 'f5'; names(fit_f5_prog$strat) # must match legend.labs
plot_f5_ <- get_plot_km(fit=fit_f5_prog, data=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), title='               ', ylab='PFS\nprobability',
                          leg.labs=c('enhancing < 1/3', 'enhancing >= 1/3'),  pal=c(wes_palette("Rushmore1")[3], wes_palette("Royal1")[1])); print(plot_f5_)

cn <- 'f6'; names(fit_f6_prog$strat) 
plot_f6_ <- get_plot_km(fit=fit_f6_prog, data=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), title='               ', ylab='PFS\nprobability',
                          leg.labs=c('nCET < 1/3', 'nCET >= 1/3'),  pal=c(wes_palette("Rushmore1")[3], wes_palette("Royal1")[1])); print(plot_f6_)

cn <- 'f7'; names(fit_f7_prog$strat) 
plot_f7_ <- get_plot_km(fit=fit_f7_prog, data=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), title='               ', ylab='PFS\nprobability',
                          leg.labs=c('necrosis < 1/3', 'necrosis >= 1/3'),  pal=c(wes_palette("Rushmore1")[3], wes_palette("Royal1")[1])); print(plot_f7_)

cn <- 'f9'; names(fit_f9_prog$strat) 
plot_f9_ <- get_plot_km(fit=fit_f9_prog, data=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), title='               ', ylab='PFS\nprobability',
                          leg.labs=c('focal', 'non-focal'),  pal=c(wes_palette("Rushmore1")[3], wes_palette("Royal1")[1])); print(plot_f9_)

cn <- 'f10'; names(fit_f10_prog$strata) 
plot_f10_ <- get_plot_km(fit=fit_f10_prog, data=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), title='               ', ylab='PFS\nprobability',
                           leg.labs=c('expansive', 'infiltrative'), pal=c(wes_palette("Rushmore1")[3], wes_palette("Royal1")[1])); print(plot_f10_)

cn <- 'f14'; names(fit_f14_prog$strat) 
plot_f14_ <- get_plot_km(fit=fit_f14_prog, data=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), title='               ', ylab='PFS\nprobability',
                           leg.labs=c('edema < 1/3', 'edema >= 1/3'),  pal=c(wes_palette("Rushmore1")[3], wes_palette("Royal1")[1])); print(plot_f14_)
a <- plot_grid(plot_f5_$plot, plot_f5_$table, nrow=2, rel_heights=c(1, .4))
b <- plot_grid(plot_f6_$plot, plot_f6_$table, nrow=2, rel_heights=c(1, .4))
c <- plot_grid(plot_f7_$plot, plot_f7_$table, nrow=2, rel_heights=c(1, .4))
d <- plot_grid(plot_f9_$plot, plot_f9_$table, nrow=2, rel_heights=c(1, .4))
e <- plot_grid(plot_f10_$plot, plot_f10_$table, nrow=2, rel_heights=c(1, .4))
f <- plot_grid(plot_f14_$plot, plot_f14_$table, nrow=2, rel_heights=c(1, .4))

img_prog <- plot_grid(a, b, c, d, e, f, ncol=1)

plot_grid(img_os, NULL, img_prog, nrow=1, rel_widths=c(1, 0.1, 1), labels=c('a','','b'), label_size=16)
ff <- file.path(fdd, 'km_img_individual.pdf'); print(ff)
dev.print(pdf, ff, height=20, width=8)


#
# COMBINE radiogenomic kms ---------

# f6 and chr3p21 surv
cn <- 'f6'; df <- chromosome[[cn]]; colnames(df); names(fit_f6_surv_chr3p21$strata)
plot_f6_chr <- get_plot_km(fit=fit_f6_surv_chr3p21, data= left_join(df[df$X %in% vaskm$X, c('X','chr3p21')], vaskm), # SURVIVAL DATA
                           title='               ', leg.labs=c('nCET + chr3p21 not enriched',
                                                               'nCET + chr3p21 enriched'), 
                           pal=wes_palette("BottleRocket2")[c(1,2)]); print(plot_f6_chr)

# f10 surv
cn <- 'f10'; df <- chromosome[[cn]]; colnames(df) ; names(fit_f10_chr19p13$strata)
plot_f10_chr19 <- get_plot_km(fit=fit_f10_chr19p13, data=left_join(df[df$X %in% vaskm$X, c('X','chr19p13')], vaskm), # SURVIVAL DATA
                              title='               ', leg.labs=c('infiltrative + chr19p13 not enriched',
                                                                  'infiltrative + chr19p13 enriched'),
                              pal=wes_palette("BottleRocket2")[c(1,2)]); print(plot_f10_chr19)
names(fit_f10_chr1p35$strata)
plot_f10_chr1 <- get_plot_km(fit=fit_f10_chr1p35, data=left_join(df[df$X %in% vaskm$X, c('X','chr1p35')], vaskm), # SURVIVAL DATA
                             title='               ', leg.labs=c('infiltrative + chr1p35 not enriched',
                                                                 'infiltrative + chr1p35 enriched'),
                             pal=wes_palette("BottleRocket2")[c(1,2)]); print(plot_f10_chr1)

# replot os for vasari
plot_os_vasari <- get_plot_km(fit=fit_vasari_os, data=vaskm, # SURVIVAL DATA
                              title='               ', leg.labs=c('VASARI cohort'),
                              pal=wes_palette("Rushmore1")[3]); print(plot_os_vasari)

# combine surv plots
# a <- plot_grid(plot_os_vasari$plot, plot_os_vasari$table, nrow=2, rel_heights=c(1, .3))
b <- plot_grid(plot_f6_os$plot, plot_f6_os$table, nrow=2, rel_heights=c(1, .32))
c <- plot_grid(plot_f6_chr$plot, plot_f6_chr$table, nrow=2, rel_heights=c(1, .32))
d <- plot_grid(plot_f10_os$plot, plot_f10_os$table, nrow=2, rel_heights=c(1, .32))
e <- plot_grid(plot_f10_chr19$plot, plot_f10_chr19$table, nrow=2, rel_heights=c(1, .32))
f <- plot_grid(plot_f10_chr1$plot, plot_f10_chr1$table, nrow=2, rel_heights=c(1, .32))

partf6 <-  plot_grid(b, NULL, c, ncol=3, rel_widths=c(1, 0.05, 1))
partf10 <- plot_grid(d, NULL, e, ncol=3, rel_widths=c(1, 0.05, 1))
partf10_b <- plot_grid(NULL, NULL, f, ncol=3, rel_widths=c(1, 0.05, 1))

plot_os_comb <- plot_grid(partf6, NULL, partf10, NULL, partf10_b, ncol=1, rel_heights=c(1, 0.03, 1, 0.03, 1))
plot_os_comb
#rel_heights=c(1,0.1,1,1), labels=c('a','','b'), label_size=16, vjust=0.9)
# progression use ggsurvplot----------

# just f6 
cn <- 'f6'; names(fit_f6_prog$strata) # print fit to get median survival
plot_f6_prog <- get_plot_km(fit=fit_f6_prog, data=bind_cols(labels[labels$X %in% vasprog$X, c('X',cn)], vasprog), # PROGRESSION DATA
                            title='               ', leg.labs=c('nCET < 1/3', 
                                                                'nCET >= 1/3'),  ylab='PFS\nprobability',
                            pal=c(wes_palette("Royal1")[1], wes_palette("Rushmore1")[3])); print(plot_f6_prog)

# f6 and nl curves
cn <- 'f6'; df <- verhaak_counts[[cn]]; colnames(df); names(fit_f6_prog_NL$strat)
plot_f6_nl <- get_plot_km(fit=fit_f6_prog_NL, data= left_join(df[df$X %in% vasprog$X, c('X','NL')], vasprog), # PROGRESSION DATA
                          title='               ', leg.labs=c('nCET + NL not enriched',
                                                              'nCET + NL enriched'),
                          pal=c(wes_palette("Rushmore1")[3], wes_palette("BottleRocket2")[1]), ylab='PFS\nprobability'); print(plot_f6_nl)

# combine prog plots
a <- plot_grid(plot_f6_prog$plot, plot_f6_prog$table, nrow=2, rel_heights=c(1, .32))
b <- plot_grid(plot_f6_nl$plot, plot_f6_nl$table, nrow=2, rel_heights=c(1, .32))

get_km_prog_comb <- plot_grid(a, NULL, b, ncol=3, rel_widths=c(1, 0.05, 1), labels=c('a','','b'), label_size=16, vjust=1)
get_km_prog_comb

#
# combine surv and progression use ggsurvplot----------
plot_grid(get_km_prog_comb, NULL , plot_os_comb, nrow=3, rel_heights=c(0.33, 0.015, 1))
ff <- file.path(fdd, 'sig_radiogenomic_km.pdf'); print(ff)
dev.print(pdf, ff, height=15, width=8.75)

#
# replot TCGA-GBM surv and progression use ggsurvplot----------
# calc values
summary(fit_all_os)$table

# TCGA-GBM 
plot_all_os <- get_plot_km(fit=fit_all_os, data=km,  title='               ', leg.labs=NA, 
                           pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_all_os)
# vasari
plot_vas_os <- get_plot_km(fit=fit_vasari_os, data=vaskm,  title='               ', leg.labs=NA, 
                           pal=c(wes_palette("Royal1")[1], wes_palette("BottleRocket2")[2])); print(plot_vas_os)

# radiognomic
rg <- dfs
rg$cohort <- 'radiogenomic'
rg <- rg[, c('X','cohort')]
rg <- left_join(rg, vaskm[, c('X','outcome','yrs')])
survfit(Surv(rg$yrs, rg$outcome) ~ 1);

# vasari vs TCGA-GBM surv
vs <- vaskm[, c('X','outcome','yrs')]
vs$cohort <- 'Vasari'
temp <- km[, c('X','outcome','yrs')]
temp$cohort <- 'TCGA-GBM'

vs <- bind_rows(vs[, colnames(temp)], temp)
vs <- bind_rows(rg[, colnames(vs)], vs)
vs$cohort <- factor(vs$cohort, levels=c('TCGA-GBM', 'Vasari','radiogenomic'))


cn <- 'cohort'; fit_vs_os <- get_km_s(df=vs, y_colname=cn, cn); fit_vs_os
plot_vs_os <- get_plot_km(fit=fit_vs_os, data=vs,  title='               ', leg.labs=c('TCGA-GBM','VASARI subset', 'radiogenomic'), 
                          xlim=c(0,11), ylab='OS\nprobability', 
                          pal=c(wes_palette("Zissou1")[c(2,3)], wes_palette("BottleRocket2")[2])); print(plot_vs_os)

# vasari vs TCGA-GBM prog
vs <- vasprog[, c('X','outcome','yrs')]
vs$cohort <- 'Vasari'
temp <- prog[, c('X','outcome','yrs')]
temp$cohort <- 'TCGA-GBM'
rg <- dfs
rg$cohort <- 'radiogenomic'
rg <- rg[, c('X','cohort')]
rg <- left_join(rg, vasprog[, c('X','outcome','yrs')])
survfit(Surv(rg$yrs, rg$outcome) ~ 1);

vs <- bind_rows(vs[, colnames(vs) %in% colnames(temp)], temp)
vs <- bind_rows(rg[, colnames(vs)], vs)
vs$cohort <- factor(vs$cohort, levels=c('TCGA-GBM', 'Vasari','radiogenomic'))

cn <- 'cohort'; fit_vs_prog <- get_km_p(df=vs, y_colname=cn, cn) ; fit_vs_prog
plot_vs_prog <- get_plot_km(fit=fit_vs_prog, data=vs,  title='               ', leg.labs=c('TCGA-GBM','VASARI subset','radiogenomic'), 
                            xlim=c(0,11), ylab='PFS\nprobability', 
                            pal=c(wes_palette("Zissou1")[c(2,3)], wes_palette("Rushmore1")[3])); print(plot_vs_prog)

# combine vasari vs TCGA-GBM plots
a <- plot_grid(plot_vs_os$plot, plot_vs_os$table, nrow=2, rel_heights=c(1, .32))
b <- plot_grid(plot_vs_prog$plot, plot_vs_prog$table, nrow=2, rel_heights=c(1, .32))

plot_grid(a, NULL, b, ncol=3, rel_widths=c(1, 0.05, 1), labels=c('a','b'), label_size=16, vjust=1)
ff <- file.path(fdd, 'km_all_vs_vasari.pdf'); print(ff)
dev.print(pdf, ff, height=4, width=8.5)


#


#
# took at patient timelines with diff km curves ----------

# clean drug info
events$type[events$type %in% c('[not available]', 'other')] <- '-1'
events$type[events$type %in% c('gliadel bcnu', 'gliadel wafer', 'gliadel wafers', 'gliadel waters', 'gliadel', 
                               'carmustine (bcnu)')] <- 'carmustine'
events$type[events$type %in% c('dexamethasome', 'dexamethazone', 'dexamethsone', 'dexmethasone', 'dexaethasone')] <- 'dexamethasone'
events$type[events$type %in% c('hyroxyurea', 'hydroxurea')] <- 'hydroxyurea'
events$type[events$type %in% c('cpt 11', 'cpt-11')] <- 'irinotecan'
events$type[events$type %in% c('tamoxiten')] <- 'tamoxifen'
events$type[events$type %in% c('cisplatain', 'cddp')] <- 'cisplatin'
events$type[events$type %in% c('cilenaitide')] <- 'cilengitide'
events$type[events$type %in% c('temador', 'temodar', 'temodor', 'temoxolomide', 'temozolamide', 'temozolomoide',
                               'temozolomode', 'temozolomide', 'temozomide',
                               'metronomic temodar')] <- 'temozolomide'
events$type[events$type %in% c('tarceva', 'tanceva')] <- 'erlotinib'
events$type[events$type %in% c('cai (nabtt 9712)', 'cai (nabtt 97212)', 'cai nabit 9712')] <- 'cai'
events$type[events$type %in% c('avastin', 'bevacizumab or placebo rtog 0825')] <- 'bevacizumab'
events$type[events$type %in% c('arsenic trioxide (ato)', 'arsenic tnoxide')] <- 'arsenic trioxide'
events$type[events$type %in% c('tarceva', 'carmustine (bcnu)', 'bcnu')] <- 'carmustine'
events$type[events$type %in% c('cra')] <- 'cis-retinoic acid'
events$type[events$type %in% c('tipifarnib (r115777)', 'tipfarnib (r115777)')] <- 'tipifarnib'
events$type[events$type %in% c('ci980', 'ci 980')] <- 'ci-980'
events$type[events$type %in% c('mab i-131', 'mab i131', 'mab i 131')] <- 'mabi131'
events$type[events$type %in% c('bsi-201', 'bs1-201')] <- 'iniparib'
events$type[events$type %in% c('levenracetam')] <- 'levetiracetam'
events$type[events$type %in% c('lumustine', 'ccnu')] <- 'lomustine'
events$type[events$type %in% c('motexatin gadoinium', 'metexafin gadolinium')] <- 'motexafin gadolinium'
events$type[events$type %in% c('vp-16')] <- 'etoposide'
events$type[events$type %in% c('paclitaxel')] <- 'taxol'

sort(unique(events$type))
length(unique(events$type))

z <- table(events$type)
z <- as.data.frame(z)

# progression
tm <- full_join(km, events)
