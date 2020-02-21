#' Analyze the bootstrapped model's results
#' 
library(dplyr)
library(reshape2)
ddir <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/experiments'
fdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_figs/bootstrap'

# if previously saved, load observed classification performance plots
bin_plot <- readRDS(file=file.path(getwd(), 'observed_cv_performance_binary.rds'))
mul_plot <-  readRDS(file=file.path(getwd(), 'observed_cv_performance_multi.rds'))


# bootstrap method: do cv split --> for each split --> resample and remove samples that do not have at least 1 example of each class in val fold
dd1 <- file.path(ddir,'boot_others_19_01_22_strat_cv_new_seeds')
dd1_nn <- file.path(ddir,'gbm_ae_19_04_15_vasari_boot')


# read boot data from other models or nn models
get_data <- function(dir_path, models, labels, boots, others, pred_type){
  # dir_path <- dd1
  # # models <- c('svmc', 'gbc', 'logit2', 'rfc', 'logit1')
  # model_name <- 'logit2'
  # label_name <- 'f10'
  bs <- lapply(models, function(model_name) {cat(model_name)
    t <- lapply(labels, function(label_name) { cat(label_name)
      a <- lapply(boots, function(i) { #cat(i)
        if (others){
          # read other model's boots
          res <- file.path(dir_path, i, model_name, label_name, 'cv_results.csv')
          if (file.exists(res)) read.csv(res, header=T)$mean_test_score
        } else {
          # read nn boots
          if (pred_type=='binary') {
            fn = 'cv_logger.txt'
            metric = 'val_roc'
          } else {
            fn = 'cv_scores.txt'
            metric = 'val_f1_micro'
          }
          res <- list.files(path=file.path(dir_path, i, label_name, model_name), pattern=fn, recursive=T, full.names=T)
          if (file.exists(res)) read.csv(res, header=T, sep='\t', row.names=1)['avg',metric]
        }
        })
      unlist(a)
    })
    # if no cv results for a label, remove them
    r <- unlist(lapply(t, function(i) !is.null(i))) 
    t <- t[r]
    
    if (length(t)==0) { # if no cv results for any label
      print(length(t))
      t <- NULL
    } else {
      names(t) <- labels[r]
      # format
      t <- as.data.frame(t)
      t$model <- model_name
      row.names(t) <- boots
      t$boots <- boots
      t <- melt(t)
    }
    t
  })
}

# read all metric data from boots for nn models
get_data_nn_all <- function(dir_path, labels, boots, pred_type){
  # dir_path <- dd1_nn
  # label_name <- 'f10'
  # 
  model_name <- 'neuralnets'
  t <- lapply(labels, function(label_name) { cat(label_name)
    a <- lapply(b, function(i) { #cat(i)
      # read nn boots
      res <- list.files(path=file.path(dir_path, i, label_name, model_name), pattern='cv_logger.txt', recursive=T, full.names=T)
      res <- if (file.exists(res)) read.csv(res, header=T, sep='\t', row.names=1)
      res$fold <- rownames(res)
      res$boot <- i
      res$label <- label_name
      if (pred_type=='multiClass') {
        f1 <- list.files(path=file.path(dir_path, i, label_name, model_name), pattern='cv_scores.txt', recursive=T, full.names=T)
        f1 <- if (file.exists(f1)) read.csv(f1, header=T, sep='\t', row.names=1)[, c('train_f1_micro','val_f1_micro')]
        f1$fold <- rownames(f1)
        res <- left_join(res, f1, by='fold')
      } 
      res
    })
    a <- bind_rows(a)
  })
}

# get the distribution of individual model's boot scores
get_dist <- function(dd){
  boot_aucs <- get_data(dir_path=dd,  models=models, labels=labels, boots=boots)
  boot_aucs <- bind_rows(boot_aucs)
  colnames(boot_aucs) <- c('label','auc','model')
  
  
  # ints and slopes for qq plot
  # see, https://mgimond.github.io/ES218/Week06a.html
  intsl <- boot_aucs %>% group_by(label, model) %>% 
    summarize(q25    = quantile(auc, 0.25, type=5),
              q75    = quantile(auc, 0.75, type=5),
              norm25 = qnorm(0.25),
              norm75 = qnorm(0.75),
              slope  = (q25 - q75) / (norm25 - norm75),
              int    = q25 - slope * norm25) %>%
    select(model, label, slope, int)
  
  
  # norm test pvalues
  boot_pvalues <- boot_aucs %>% group_by(label, model) %>% 
    summarize( p.value = shapiro.test(auc)$p.value ) %>%
    select(model, label, p.value)
  
  boot_pvalues$plabel <- paste0('p = ', round(boot_pvalues$p.value, digits=3))
  boot_pvalues$x <- 1
  boot_pvalues$y <- 0.3
  
  
  # plot qq
  # ggplot(boot_aucs, aes(sample=auc)) +
  #   stat_qq(shape=1) +
  #   # geom_abline(data=intsl, aes(intercept=int, slope=slope), col='red') +
  #   geom_text(data=boot_pvalues, aes(x=x, y=y, label=plabel), inherit.aes=F) +
  #   theme_bw() +
  #   facet_grid(label~model)
  
  # plot hist
  boot_pvalues$x <- .2
  boot_pvalues$y <- 60
  ggplot(boot_aucs, aes(x=auc)) +
    geom_histogram(bins=50)+
    # geom_abline(data=intsl, aes(intercept=int, slope=slope), col='red') +
    geom_text(data=boot_pvalues, aes(x=x, y=y, label=plabel), inherit.aes=F) +
    theme_bw() +
    facet_grid(label~model)
  
  # confidence intervials
  boot_cis <- boot_aucs %>% group_by(label, model) %>% 
    summarize(q025 = quantile(auc, 0.025, type=5),
              q975 = quantile(auc, 0.975, type=5)) %>%
    select(model, label, q025, q975)
  boot_cis <- as.data.frame(boot_cis)
  
  boot_cis2 <- boot_aucs %>% group_by(label, model) %>% 
    summarize(q025 = quantile(auc, 0.025, type=5),
              q975 = quantile(auc, 0.975, type=5)) %>%
    select(model, label, q025, q975)
  boot_cis <- as.data.frame(boot_cis)
  
  # write data out
  write.csv(boot_cis, file.path(dd, 'conf_intervals.csv'), row.names=F)
  return(boot_cis)
}


# ci_1 <- get_dist(dd1)
# ci_2 <- get_dist(dd2)

models <-  c('svmc', 'gbc', 'logit2', 'rfc', 'logit1') # note that logit1 is selected for f7, and logit2 for all others
labels <- c('f5','f6','f7','f9','f10','f14')
b <- paste0('boot_', seq(0,99)) # boots to check

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# invidiual  ---------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

# read in boots
boots_nn <- get_data(dir_path=dd1_nn, models=c('neuralnets'), labels=labels, others=F, pred_type='binary', boots=b)[[1]] # binary

boots_ot <- get_data(dir_path=dd1, models=c('svmc', 'gbc', 'rfc'), labels=labels, others=T, pred_type='binary', boots=b) # binary
boots_ot_b1 <- get_data(dir_path=dd1, models=c('logit2'), labels=c('f5','f6','f9','f10','f14'), others=T, pred_type='binary', boots=b) # binary
boots_ot_b2 <- get_data(dir_path=dd1, models=c('logit1'), labels=c('f7'), others=T, pred_type='binary', boots=b) # binary
boots_ot <- bind_rows(boots_ot, boots_ot_b1, boots_ot_b2) 

boots_ot$model[boots_ot$model %in% c('logit1','logit2')] <- 'logit' # rename for plotting
boots_ot$model[boots_ot$model %in% c('gbc')] <- 'gbt'  
boots_ot$model[boots_ot$model %in% c('rfc')] <- 'rf' 
boots_ot$model[boots_ot$model %in% c('svmc')] <- 'svm' 
boots_nn$model[boots_nn$model %in% c('neuralnets')] <- 'nn'

rm(boots_ot_b1, boots_ot_b2)

# individual model metrics:
boots_ind <- bind_rows(boots_ot, boots_nn)
colnames(boots_ind) <- c('model', 'boots', 'label','score')

boots_ind$model <- factor(boots_ind$model, levels=rev(c('logit','svm','rf','gbt','nn'))) # refactor for plotting
boots_ind$label <- mapvalues(boots_ind$label, from = c('f5', 'f6','f7','f14','f10','f9'), 
                             to = c('enhancing',   'nCET', 'necrosis', 'edema',      'infilatrative',    'focal'))
boots_ind$label <- factor(boots_ind$label, levels=c('enhancing',   'nCET', 'necrosis', 'edema',      'infilatrative',    'focal'))


# 95% CI
boots_ind_cis <- boots_ind %>% group_by(model, label) %>%
  dplyr::summarize(q025 = round(quantile(score, 0.025, type=5), digits=3),
            q975 = round(quantile(score, 0.975, type=5), digits=3)) %>%
  dplyr::select(model, label, q025, q975)

boots_ind_cis$text <- paste0('(', boots_ind_cis$q025, ',', boots_ind_cis$q975,')')
boots_ind_cis$x <- 0.8
boots_ind_cis$y <- 15

# norm test pvalues
boot_pvalues_ind <- boots_ind %>% group_by(model, label) %>% 
  dplyr::summarize( p.value = shapiro.test(score)$p.value ) %>% # P < 0.5 = IS NOT NORMAL
  dplyr::select(model, label, p.value)
boot_pvalues_ind$plabel <- paste0('p = ', round(boot_pvalues_ind$p.value, digits=3))
boot_pvalues_ind$x <- .8
boot_pvalues_ind$y <- 15

saveRDS(boots_ind, file.path(getwd(), 'plot_data', 'boots_ind.rds'))
saveRDS(boot_pvalues_ind, file.path(getwd(), 'plot_data', 'boot_pvalues_ind.rds'))
saveRDS(boots_ind_cis, file.path(getwd(), 'plot_data', 'boots_ind_cis.rds'))


#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# difference: nn - other models  ---------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

# get the mean difference between nn and other models boot scores
boots <- lapply(unique(boots_ot$model), function(i) { cat(i)
  a <- lapply(unique(boots_nn$variable), function(l) { cat(l)
    nn <- boots_nn[boots_nn$variable==l,]
    ot <- boots_ot[(boots_ot$variable==l &  boots_ot$model==i),]
    if ( all.equal(nn$boots, ot$boots) ) { # check boots are the same rows
      d <- nn[, c('variable'), drop=F]
      d[, 'difference'] <- nn$value - ot$value
      d[, -1, drop=F]
      d$model <- i
      d$boots <- nn$boots
      d
    }
  })
  a <- bind_rows(a)
})
boots <- bind_rows(boots)
colnames(boots) <- c('label','difference','model','boots')
boots$model <- factor(boots$model, levels=rev(c('logit','svm','rf','gbt','nn'))) # refactor for plotting

boots$label <- mapvalues(boots$label, from = c('f5', 'f6','f7','f14','f10','f9'), 
                      to = c('enhancing', 
                             'nCET',
                             'necrosis',
                             'edema', 
                             'infilatrative',
                             'focal'))
boots$label <- factor(boots$label, levels=c('enhancing', 
                                            'nCET',
                                            'necrosis',
                                            'edema', 
                                            'infilatrative',
                                            'focal'))

# boots_qq <- boots %>% group_by(model, label) %>% 
#   summarize(q25    = quantile(value, 0.25, type=5),
#             q75    = quantile(value, 0.75, type=5),
#             norm25 = qnorm(0.25),
#             norm75 = qnorm(0.75),
#             slope  = (q25 - q75) / (norm25 - norm75),
#             int    = q25 - slope * norm25) %>%
#   select(model, label, slope, int)

boot_cis <- boots %>% dplyr::group_by(model, label) %>% 
  dplyr::summarize(q025 = round(quantile(difference, 0.025, type=5), digits=3),
            q975 = round(quantile(difference, 0.975, type=5), digits=3)) %>%
  dplyr::select(model, label, q025, q975)
boot_cis <- as.data.frame(boot_cis)
boot_cis$text <- paste0('(', boot_cis$q025, ',\n', boot_cis$q975,')')
boot_cis$x <- 0.45
boot_cis$y <- 10
# norm test pvalues
boot_pvalues <- boots %>% dplyr::group_by(model, label) %>% 
  dplyr::summarize( p.value = shapiro.test(difference)$p.value ) %>% # P < 0.5 = IS NOT NORMAL
  dplyr::select(model, label, p.value)

boot_pvalues$plabel <- paste0('p = ', round(boot_pvalues$p.value, digits=3))

# save data 
saveRDS(boots, file.path(getwd(), 'plot_data', 'boots.rds'))
saveRDS(boot_pvalues, file.path(getwd(), 'plot_data', 'boot_pvalues.rds'))
saveRDS(boot_cis, file.path(getwd(), 'plot_data', 'boot_cis.rds'))



#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# nn data  ---------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------


# get nn data for all metrics values
t <- get_data_nn_all(dir_path=dd1_nn, labels=labels, boots=b, pred_type='binary')
t <- bind_rows(t)
t$label <- mapvalues(t$label, from = c('f5', 'f6','f7','f14','f10','f9'), 
                         to = c('enhancing', 
                                'nCET',
                                'necrosis',
                                'edema', 
                                'infilatrative',
                                'focal'))
t$label <- factor(t$label, levels=c('enhancing', 
                                            'nCET',
                                            'necrosis',
                                            'edema', 
                                            'infilatrative',
                                            'focal'))
# save data 
saveRDS(t, file.path(getwd(), 'plot_data', 'boots_nn_all.rds'))
