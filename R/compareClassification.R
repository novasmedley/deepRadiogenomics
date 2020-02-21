#---------------------------------------- Description---------------------------------------- 
#' Compare the classification performance between other models and nn

#'---------------------------------------- ---------------------------------------- 
library(dplyr)
library(tidyr)
library(reshape2) # melt
library(ggplot2)

dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/experiments'
fdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/temp_figs'

# autoencoder param search  ------
ae <- read.csv(file.path(dd, 'gbm_ae_19_04_15_cv', 'parsed_cv_results.csv'), stringsAsFactors=F)

ae$opt_act <- paste0(ae$opt, '_',ae$act)
ggplot(ae, aes(x=r2, y=val_r2, color=archit, shape=opt_act)) +
  geom_point(size=3) + theme_bw() +theme(axis.title.y = element_text(angle=0, vjust=0.5))
ff <- file.path(fdd, 'ae_compare_hyperparams.png');print(ff)
dev.print(png, ff, res=300, height=4, width=6, units="in")

# autoencoder retrain metrics ------
ae <- read.csv(file.path(dd, 'gbm_ae_19_04_15_retrain', 'autoencoder','neuralnets',
                         '4000_2000_1000_0_0_tanh_decay_0_drop_0_opt_Adadelta_loss_mae_bat_50_eph_500',
                         'retrain_metrics.txt'), stringsAsFactors=F, sep='\t')

ggplot(ae, aes(x=r2)) +
  geom_histogram(color='black',fill='white') +
  theme_grey() + theme(axis.title.y = element_text(angle=0, vjust=0.5))
ff <- file.path(fdd, 'ae_retrained_r2_hist.png');print(ff)
dev.print(png, ff, res=300, height=4, width=6, units="in")



# compare cv performance ae based models  ------
exp <- 'gbm_ae_19_04_15_vasari_cv'

ae_based_nn <- read.csv(file.path(dd, exp, 'parsed_cv_results_binary.csv'))
ae_based_nn$type <- unlist(lapply(ae_based_nn$nl, function(i) ifelse(i=='[]', 'vanilla','pretrained')))
a <- ae_based_nn %>% group_by(label,type) %>%
  top_n(1, val_roc) 


ggplot(a, aes(x=type, y=val_roc, fill=type)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black", lwd=0.3)+
  theme_bw() + theme(axis.title.y = element_text(angle=0, vjust=0.5),
                     axis.text.x=element_blank()) +
  scale_y_continuous(minor_breaks=seq(.55, 1, .05), breaks=seq(.5, 1, .1), limits=c(.48,1), oob=rescale_none, expand = c(0.001,0.001)) +
  facet_grid(.~label)

ff <- file.path(fdd, 'nn_nnae_performance.png');print(ff)
dev.print(png, ff, res=300, height=3, width=7, units="in")


ggplot(a, aes(x=type, y=roc, fill=type)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black", lwd=0.3)+
  theme_bw() + theme(axis.title.y = element_text(angle=0, vjust=0.5),
                     axis.text.x=element_blank()) +
  scale_y_continuous(minor_breaks=seq(.55, 1, .05), breaks=seq(.5, 1, .1), limits=c(.48,1), oob=rescale_none, expand = c(0.001,0.001)) +
  facet_grid(.~label)

ff <- file.path(fdd, 'nn_nnae_performance_train.png');print(ff)
dev.print(png, ff, res=300, height=3, width=7, units="in")

# compare wts between ae based models  ------
exp <- 'gbm_ae_19_04_15_vasari_retrain'

f5 <- read.csv(file.path(dd, exp, 'f5_nn_nnae_correlations.csv'))
f6 <- read.csv(file.path(dd, exp, 'f6_nn_nnae_correlations.csv'))
f7 <- read.csv(file.path(dd, exp, 'f7_nn_nnae_correlations.csv'))
f5$label <- 'f5'
f6$label <- 'f6'
f7$label <- 'f7'

corrs <- bind_rows(f5,f6,f7)
rm(f5,f6,f7)

corrs <- gather(corrs, layer, corr, 2:4)

breaks <- seq(-.2, .2, 0.05)
ggplot(corrs, aes(x=X, y=corr, color=corr)) +
  geom_point(size=.5, alpha=.7) +
  xlab('node') +
  scale_color_gradientn(colours=c(viridis(5, direction=1, option = "C"), viridis(5, direction=-1, option = "C")), breaks=breaks, limits=c(-.2,.2))+
  theme_bw() + theme(axis.title.y = element_text(angle=0, vjust=0.5)) +
  facet_grid(label ~ layer, scales='free_x', space='free_x') 

ff <- file.path(fdd, 'nn_nnae_correlations.png');print(ff)
dev.print(png, ff, res=300, height=4, width=8, units="in")

# compare neural nets   ------
nn_no_pre <- 'gbm_vasari_cv_200patience_18_12_20'
nn_pre <- 'gbm_vasari_cv_200patience_18_12_20_with_ae_pretrain_grid_search_freeze'
nn_pre_ae_based <- 'gbm_vasari_aebased_19_04_11_cv'

bin_labs <- c('f5','f6','f7','f9','f10','f14')

a <- read.csv(file.path(dd, nn_no_pre, 'parsed_cv_results_binary.csv'), stringsAsFactors=F)
b <- read.csv(file.path(dd, nn_pre, 'parsed_cv_results_binary.csv'), stringsAsFactors=F)
a$model <- 'nn_no_pre'; b$model <- 'nn_with_pre'
bin_nn <- bind_rows(a,b)
bin_nn <- bin_nn[bin_nn$label %in% bin_labs, ] # get just the binary labels
unique(bin_nn$label)

a <- read.csv(file.path(dd, nn_no_pre, 'parsed_cv_results_multiclass.csv'), stringsAsFactors=F)
b <- read.csv(file.path(dd, nn_pre, 'parsed_cv_results_multiclass.csv'), stringsAsFactors=F)
a$model <- 'nn_no_pre'; b$model <- 'nn_with_pre'
mul_nn <- bind_rows(a,b)

mul_nn$label[mul_nn$label=='f10_2'] <- 'f10' # rename
mul_nn <- mul_nn[mul_nn$label=='f10', ] # # get just the multiclass label, select the f10 with f1_scores
unique(mul_nn$label) # check 
unique(mul_nn$model) # check 
 
# compare between hyperparams with nn_pre
plot_bin <- bin_nn
plot_bin$opt_act <- paste0(plot_bin$opt, '_',plot_bin$act)
ggplot(plot_bin, aes(x=roc, y=val_roc, color=drop, shape=opt_act)) +
  geom_point(size=3) + theme_bw()+ theme(axis.title.y = element_text(angle=0, vjust=0.5)) +
  facet_grid(label ~ model) 
ff <- file.path(fdd, 'nn_compare_hyperparams_actopt.png');print(ff)
dev.print(png, ff, res=300, height=8, width=10, units="in")


ggplot(plot_bin, aes(x=roc, y=val_roc, color=drop, shape=nl)) +
  geom_point(size=3) + theme_bw() + theme(axis.title.y = element_text(angle=0, vjust=0.5)) +
  facet_grid(label ~ model) 
ff <- file.path(fdd, 'nn_compare_hyperparams_freezenl.png');print(ff)
dev.print(png, ff, res=300, height=8, width=10, units="in")

ggplot(plot_bin, aes(x=roc, y=val_roc, color=drop, shape=archit)) +
  geom_point(size=3) + theme_bw() + theme(axis.title.y = element_text(angle=0, vjust=0.5)) +
  facet_grid(label ~ model) 
ff <- file.path(fdd, 'nn_compare_hyperparams_archi.png');print(ff)
dev.print(png, ff, res=300, height=8, width=10, units="in")

# compare between hyperparams with nn_pre_ae_based
a <- read.csv(file.path(dd, nn_no_pre, 'parsed_cv_results_binary.csv'), stringsAsFactors=F)
b <- read.csv(file.path(dd, nn_pre_ae_based, 'parsed_cv_results_binary.csv'), stringsAsFactors=F)
a$model <- 'nn_no_pre'; b$model <- 'nn_with_pre'
bin_nn <- bind_rows(a,b)
bin_nn <- bin_nn[bin_nn$label %in% bin_labs, ] # get just the binary labels
unique(bin_nn$label)
plot_bin <- bin_nn 

ggplot(plot_bin[plot_bin$archit=='1000_500_250_0_0' & plot_bin$opt=='Nadam' & plot_bin$act=='tanh' & plot_bin$drop<.8 ,], 
       aes(x=roc, y=val_roc, color=drop, shape=nl)) + 
  geom_point(size=3) + theme_bw() + theme(axis.title.y = element_text(angle=0, vjust=0.5)) +
  facet_grid(label ~ model) 
ff <- file.path(fdd, 'nn_compare_hyperparams_aebased_1000Nadamtanh.png');print(ff)
dev.print(png, ff, res=300, height=8, width=10, units="in")

# selection best one
bin_nn$val_err <- bin_nn$val_roc
mul_nn$val_err <- mul_nn$val_f1_micro
these <- c('label','val_err','model')
nn <- bind_rows(bin_nn[, these], mul_nn[, these])

b <- nn[, c('val_err', 'label', 'model')] # get just the best one
b <- b %>% 
  group_by(model, label) %>%
  slice(which.max(val_err))


# compare models  ------
nn_pre_ae_based <- 'gbm_ae_19_04_15_vasari_cv'
bin_labs <- c('f5','f6','f7','f9','f10','f14')
nn <- read.csv(file.path(dd, nn_pre_ae_based, 'parsed_cv_results_binary.csv'), stringsAsFactors=F)
nn$model <- 'nn'
nn <- nn[nn$nl != '[]',]
# get just the best one
nn <- nn %>% 
  dplyr::group_by(label) %>%
  dplyr::slice(which.max(val_roc))

nn <- nn[, c('val_roc', 'label', 'model')] 

# other models
others <- read.csv(file.path(dd, 'others_18_12_19', 'parsed_cv_results.csv'), stringsAsFactors=F)
# boot_cis <- read.csv(file.path(dd, 'boot_others_19_01_08_stratified', 'conf_intervals.csv'), stringsAsFactors=F)

others$val_roc <- others$mean_test_score
# a <- full_join(a, boot_cis, by=c('model','label'))
others[others$model %in% c('logit1','logit2'), 'model'] <- 'logit'


others <- others %>%  # get just the best one
  dplyr::group_by(model, label) %>%
  dplyr::slice(which.max(val_roc))


# combine all mobdels
res <- bind_rows(nn, others[, colnames(nn)])
res$model <- factor(res$model, levels=rev(c('logit','svmc','rfc','gbc','nn'))) # plot order
res$label <- factor(res$label, levels=c('f5', 'f6','f7','f14','f10','f9')) # plot order
res$name <- mapvalues(res$label, from = c('f5', 'f6','f7','f14','f10','f9'), 
                       to = c('enhancing', 
                              'nCET',
                              'necrosis',
                              'edema', 
                              'infilatrative',
                              'focal'))

# DIFFERENCES BETWEEN MODELS
diff <- spread(res, model, val_roc)
diff[, 3:ncol(diff)] <- round(diff[, 3:ncol(diff)], digits=3)
diff$diff_gbc <- diff$nn - diff$gbc
diff$diff_rfc <- diff$nn - diff$rfc
diff$diff_svmc <- diff$nn - diff$svmc
diff$diff_logit <- diff$nn - diff$logit

# save data
saveRDS(res, file=file.path(getwd(), 'plot_data','cv_results.rds'))
saveRDS(diff, file=file.path(getwd(), 'plot_data','cv_results_differences.rds'))

