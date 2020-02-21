#' plot all figures here
#' 

library(ggplot2)
library(scales) # graph scales
library(wesanderson) # graph color
library(ggrepel) # geom text reprel
library(egg)
library(cowplot) # plot_grid
library(plyr) # mapvalues
library(dplyr)

dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga'
fdd <- file.path(dd, 'with_ae_paper_figs/performance')
dev.new()

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#  constants  ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
bin_labs <- c('f5','f6','f7','f9', 'f10', 'f14')

t <- 12
t2 <- t + 0

th <- theme(legend.text=element_text(size=t2, color="grey30"),
            legend.title=element_text(size=t),
            legend.position='right',
            legend.key.size=unit(.8,'line'),
            legend.margin=margin(l = -0.15, unit='cm'),
            legend.background = element_rect(fill = "transparent", colour = "transparent"),
            axis.text.y=element_text(size=t2, color="#666666"),
            axis.text.x=element_text(size=t2,color='#666666'),
            axis.title.y = element_text(margin=margin(0,5,0,0),size=t, angle=0, vjust=0.5),
            # axis.title.x = element_text(margin=margin(5,0,0,0),size=t),
            axis.title.x = element_blank(),
            strip.text.x = element_text(color='white',size=t2,margin = margin(.05,0,.05,0, "cm")),
            strip.background =  element_rect(fill='black'),
            panel.border = element_rect(colour = "black", fill=NA, size=.3),
            # panel.grid = element_line(color='grey60'),
            panel.grid.minor = element_line(colour="grey90", size=.25),
            panel.grid.major.x = element_blank()) 

x <- scale_x_discrete(expand=c(0,0))
y <- scale_y_discrete(expand=c(0,0))

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# cv  classification ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

rename_labels <- function(labels) {
  labels <- gsub('non-contrast\nenchancement', 'nCET', labels)
  labels <- gsub('infilatrative\n', 'infilatrative', labels)
  labels <- gsub('enhancement', 'enhancing', labels)
}

res <- readRDS(file.path(getwd(), 'plot_data', 'cv_results.rds'))
diff <- readRDS(file=file.path(getwd(), 'plot_data','cv_results_differences.rds'))

res$name <- rename_labels(res$name) # rename levels
diff$name <- rename_labels(diff$name)
res$name <- factor(res$name, levels=c('enhancing','nCET','necrosis', 'edema','infilatrative','focal'))
diff$name <- factor(diff$name, levels=c('enhancing','nCET','necrosis', 'edema','infilatrative','focal'))
levels(res$name) 
levels(diff$name)

res$model <- mapvalues(res$model, 
                       from = c('nn', 'gbc','logit','rfc','svmc'), 
                       to = c('nn','gbt','logit','rf','svm'))
mod_colors = c(rev(wes_palette("Darjeeling1")))
sc <- scale_fill_manual(values=mod_colors)
cv_theme <- theme_bw() + th
cv_bar <- geom_bar(stat="identity", position=position_dodge(), color = "black", lwd=0.3)

cv_plot <- ggplot(res, aes(x=name, y=val_roc, fill=model)) +  
  scale_y_continuous(minor_breaks=seq(.55, 1, .05), breaks=seq(.5, 1, .1), limits=c(.48,1), oob=rescale_none, expand = c(0.001,0.001)) +
  cv_bar + sc + cv_theme +
  ylab('mean\n AUC') + xlab('classification task');cv_plot

dev.print(png, file.path(fdd, paste0('compare_model_cv_named.','png')), res=600, height=2.25, width=7, units="in")
# setEPS(); postscript(file.path(fdd, 'compare_model_cv.eps'), height=8, width = 8.66)
# dev.print(tiff, file.path(fdd, paste0('compare_model_cv_named.','tiff')), res=600, height=2.25, width=7, units="in")

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# bootstraps invidivual models  ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

boots_ind <- readRDS(file.path(getwd(), 'plot_data', 'boots_ind.rds'))
boots_ind_cis <- readRDS(file.path(getwd(), 'plot_data', 'boots_ind_cis.rds'))

boots_ind_vline <- data.frame(name=c('enhancing', 'nCET',  'edema',  'infilatrative', 'focal'), x=c(0.5, 0.5, 0.5, 0.5, NA, 0.5))

ggplot(boots_ind, aes(x=score, fill=model)) +
  geom_histogram(bins=75)+ 
  # geom_text(data=boots_ind_cis, aes(x=x, y=y, label=text), inherit.aes=F, size=2.5) +
  geom_vline(data=boots_ind_cis, aes(xintercept=q025), colour="#BB0000", linetype="dashed") +
  geom_vline(data=boots_ind_cis, aes(xintercept=q975), colour="#BB0000", linetype="dashed") +
  geom_vline(data=boots_ind_vline, aes(xintercept=x), colour="#BB0000") + 
  sc + 
  theme_set(theme_light(base_size = 15)) + 
  theme(axis.title.y = element_text(angle = 0, vjust = .5),
        strip.text.y = element_blank()) + 
  xlab('mean AUC') +  scale_x_continuous(breaks=seq(0.3, 1.0, 0.2), limits=c(0.3,1)) +
  facet_grid(model~label)
dev.print(png, file.path(fdd, paste0('score_hist_100','.png')), res=300, height=5, width=12, units="in")

## combine
ggplot(boots_ind, aes(x=score, fill=model)) +
  # geom_histogram(bins=75, alpha = 0.8, color='black', size=.15)+
  geom_density(alpha=0.8, alhha = 0.7, size=.15) +
  # geom_text(data=boots_ind_cis, aes(x=x, y=y, label=text), inherit.aes=F, size=2.5) +
  geom_vline(data=boots_ind_cis[boots_ind_cis$model=='nn',], aes(xintercept=q025, colour=model), linetype="dashed") +
  geom_vline(data=boots_ind_cis[boots_ind_cis$model=='nn',], aes(xintercept=q975, colour=model), linetype="dashed") +
  geom_vline(data=boots_ind_vline, aes(xintercept=x), colour="#BB0000") + 
  theme_light() + theme(axis.title.y = element_text(angle = 0, vjust = .5)) +
  xlab('mean\n AUC') +  scale_x_continuous(breaks=seq(0.2, 1, 0.1)) +
  facet_grid(label~., scales='free_y')

# boxplot present version
ggplot(boots_ind, aes(x=model, y=score, fill=model)) +
  geom_boxplot() + ylab('mean\n AUC') + xlab('') + scale_y_continuous(breaks=seq(0.2, 1.0, 0.1)) + 
  sc + 
  theme_set(theme_light(base_size = 15)) + 
  theme(axis.title.y = element_text(angle = 0, hjust = .5, vjust=.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_grid(.~label)
dev.print(png, file.path(fdd, paste0('score_box_100','.png')), res=300, height=5, width=12, units="in")

# boxplot publish version
# ggplot(boots_ind, aes(x=model, y=score, fill=model)) +
#   geom_boxplot() + ylab('score') + xlab('model') + scale_y_continuous(breaks=seq(0.3, 0.9, 0.1)) + 
#   theme_light() + theme(axis.title.y = element_text(angle = 0, hjust = .5, vjust=.5)) +
#   geom_boxplot(outlier.shape=NA) + # avoid plotting outliers twice
#   geom_jitter(position=position_jitter(width=.2, height=0), fill='grey70', color='grey30', alpha=.2, size=1) +
#   facet_grid(.~label)


#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# bootstraps diff b/w models  ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

boots <- readRDS(file.path(getwd(), 'plot_data', 'boots.rds'))
boot_cis <- readRDS(file.path(getwd(), 'plot_data', 'boot_cis.rds'))
boot_pvalues <- readRDS(file.path(getwd(), 'plot_data', 'boot_pvalues.rds'))

sc <- scale_fill_manual(values=mod_colors[2:5])
# present
boot_pvalues$x <- .4
boot_pvalues$y <- 14
ggplot(boots, aes(x=difference, fill=model)) +
  geom_histogram(bins=75)+
  geom_text(data=boot_pvalues, aes(x=x, y=y, label=plabel), size=3, inherit.aes=F) +
  geom_text(data=boot_cis, aes(x=x, y=y, label=text),size=3, inherit.aes=F) +
  geom_vline(data=boot_cis, aes(xintercept=q025), colour="#BB0000", linetype="dashed") + 
  geom_vline(data=boot_cis, aes(xintercept=q975), colour="#BB0000", linetype="dashed") +
  geom_vline(xintercept=0, colour="#BB0000") + 
  sc + 
  theme_set(theme_light(base_size = 15)) + 
  theme(axis.title.y = element_text(angle = 0, vjust = .5),
        strip.text.y = element_blank()) +
  xlab('score difference\n(nn - another model)') +  #scale_x_continuous(breaks=seq(0.3, 0.9, 0.1)) +
  facet_grid(model~label) 
# dev.print(png, file.path(fdd, paste0('score_difference_hist','.png')), res=300, height=5, width=12, units="in")

ggplot(boots, aes(x=difference, fill=model)) +
  geom_histogram(bins=75, alpha = 0.8, color='black', size=.15) +
  # geom_text(data=boot_pvalues, aes(x=x, y=y, label=plabel), size=3, inherit.aes=F) +
  # geom_text(data=boot_cis, aes(x=x, y=y, label=text),size=3, inherit.aes=F) +
  # geom_vline(data=boot_cis, aes(xintercept=q025), colour="#BB0000", linetype="dashed") + 
  # geom_vline(data=boot_cis, aes(xintercept=q975), colour="#BB0000", linetype="dashed") +
  geom_vline(xintercept=0, colour="#BB0000") + 
  theme_light() + 
  theme(axis.title.y = element_text(angle = 0, vjust = .5)) +
  xlab('score difference\n(nn - another model)') +  #scale_x_continuous(breaks=seq(0.3, 0.9, 0.1)) +
  facet_grid(.~label) 

# present
ggplot(boots, aes(x=model, y=difference, fill=model)) +
  geom_boxplot() + ylab('score difference\n(nn - another model)') + xlab('label') +  theme_light() + 
  geom_hline(yintercept=0, colour="#BB0000", linetype="solid") + 
  theme(axis.title.y = element_text(angle = 0, hjust = .5, vjust=.5)) +
  scale_y_continuous(breaks=seq(-.1, .6, 0.1), limits=c(-.1,.51)) +
  facet_grid(.~label) 
dev.print(png, file.path(fdd, paste0('score_difference_box_100','.png')), res=300, height=5, width=12, units="in")

# publish
boots$label <- factor(boots$label, levels=c('enhancing', 'nCET', 'necrosis', 'edema',   'infilatrative', 'focal'))
boot_cis$label <- factor(boot_cis$label, levels=c('enhancing', 'nCET', 'necrosis', 'edema',   'infilatrative', 'focal'))

bt_theme <-  theme_bw() + th + theme(axis.title.y = element_text(angle = 0, hjust = .5, vjust=.5),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank()) 
bt_hline <-  geom_hline(yintercept=0, colour="#BB0000", linetype="solid")
boots_sc <- scale_fill_manual(values=c(wes_palette("Darjeeling2")[2], rev(wes_palette("Darjeeling1")))[3:6])
boots_sc <- scale_fill_manual(values=c("#FFEFCA",  "#EDA16A" ,"#C83741", "#6C283D", "#62BF94"),
                              labels=c('nn - gbt','nn - rf','nn - svm','nn - logit'),
                              name=' comparison')
bt_bp <- geom_boxplot(lwd=0.3, outlier.size=.65) 

boots_plot <- ggplot(boots, aes(x=model, y=difference, fill=model)) +
  geom_hline(yintercept=0.1, linetype='dashed', color='black', size=0.2) +
  geom_hline(yintercept=0, linetype='solid', color='black', size=0.2) +
  bt_bp+ boots_sc + bt_theme  + xlab('label') + 
  geom_point(data=boot_cis, aes(y=q025, color='red'), shape=8, size=.75) +
  geom_point(data=boot_cis, aes(y=q975, color='red'), shape=8, size=.75) +
  scale_color_manual(labels = c('95% CI'), values = c('red')) +
  scale_y_continuous(minor_breaks=seq(-.1,.6,0.1), breaks=seq(-.1, .6, 0.1), 
                     limits=c(-.05,.6), oob=rescale_none, expand = c(0.001,0.001)) +
  ylab('delta\nmean\nAUC') +   
  guides(fill = guide_legend(override.aes = list(shape = NA)), 
            colour=guide_legend(title='confidence\ninterval')) +
  facet_grid(.~label) ;boots_plot

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# bootstraps diff b/w models random forest ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

z <- boots_ind

z2 <- spread(z[z$model=='nn',], key=model, value=score) # nn scores
z2 <- left_join(z[z$model != 'nn',], z2) # pair nn scores with every other model's scores
z2 <- left_join(z2, boots[boots$model != 'nn', ]) # pair with differences between nn and other model's scores
z3 <- z2
z3 <- z3[z3$difference < 0, ]
z2$difference <- abs(z2$difference)

p1 <- ggplot(z2, aes(x=nn, y=score, fill=difference)) +
  geom_point(colour="black", pch=21, size=2, alpha=.7) + 
  geom_abline(slope=1, intercept=0) +
  scale_fill_viridis(limits=c(0,.6), breaks=seq(0, 0.6, 0.2),
                     option='plasma', direction=-1) +
  theme_bw() + th + theme(strip.text.y = element_text(color='white',size=t2,margin = margin(.05,0,.05,0, "cm")),
                          axis.title.x = element_text(size=t2)) +
  scale_y_continuous(breaks=seq(0.2, 1, 0.2), limits=c(0.3,1.0)) +
  scale_x_continuous(breaks=seq(0.2, 1, 0.2), limits=c(0.3,1.0)) +
  labs(y='other\nmodel\nmean\nAUC', x='neural network mean AUC', fill='difference') +
  facet_grid(label~model); p1

# cases with worse nn
# p2 <- ggplot(z3, aes(x=nn, y=score, fill=difference)) +
#   geom_point(colour="black", pch=21, size=3, alpha=.6) + 
#   geom_abline(slope=1, intercept=0) +
#   scale_x_continuous(breaks=c(seq(0.6, 0.8, 0.1))) +
#   scale_y_continuous(breaks=c(seq(0.6, 0.8, 0.1))) +
#   scale_fill_viridis(limits=c(-.1,0),option='plasma', direction=1) +
#   theme_bw() + th + theme(strip.text.y = element_text(color='white',size=t2,margin = margin(.05,0,.05,0, "cm")),
#                           axis.title.x = element_text(size=t2)) +
#   labs(y='other\nmodel\nscore', x='nn score', fill='difference') +
#   facet_grid(label~model); p2
# 
# plot_grid(p1, p2, nrow=2, rel_heights=c(1,.3), labels=c('a','b'), label_size=16)

dev.print(png, file.path(fdd, paste0('performance_100_diagonals','.png')), res=300, height=12, width=10, units="in")

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# figure performance (observed cv + bootstraps diff b/w models  ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
cv_plot2 <- plot_grid(cv_plot, NULL, nrow=1, rel_widths=c(1,.03)); cv_plot2

plot_grid(cv_plot2, boots_plot,
          labels=c('a', 'b'),
          rel_heights=c(1, 1),
          nrow=2,
          label_size=25)
dev.print(png, file.path(fdd, paste0('performance_100','.png')), res=300, height=6, width=12, units="in")



#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# bootstraps nn all metric info  ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
nn <- readRDS(file.path(getwd(), 'plot_data', 'boots_nn_all.rds'))
nn <- bind_rows(nn)
nn$label <- mapvalues(nn$label, from = c('f5', 'f6','f7','f14','f10','f9'), 
                     to = c('enhancing', 'nCET', 'necrosis', 'edema',   'expansive', 'focal'))
num_boots <- 100

t_nn <- nn[, c('label','fold','boot','epoch','acc','loss','pr','roc')]
v_nn <- nn[, c('label','fold','boot','epoch','val_acc','val_loss','val_pr','val_roc')]
colnames(t_nn) <- c('label','fold','boot','epoch','acc','loss','pr','roc')
colnames(v_nn) <- c('label','fold','boot','epoch','acc','loss','pr','roc')
t_nn$partition <- 'training'
v_nn$partition <- 'validation'

avg_t <- t_nn[t_nn$fold=='avg',]; t_nn <- t_nn[!(t_nn$fold=='avg'),]
avg_v <- v_nn[v_nn$fold=='avg',]; v_nn <- v_nn[!(v_nn$fold=='avg'),]
nn <- bind_rows(t_nn, v_nn); rm(t_nn, v_nn)

nn_theme <- theme_set(theme_light(base_size = 15))  + th + theme(strip.text = element_text(colour = 'white'),
                                    strip.text.y= element_blank(),
                                    axis.title.x = element_text(size=t),
                                    axis.text.x = element_text(size=t-2),
                                    axis.text.y = element_text(size=t-2),
                                    axis.title.y=element_text(size=t, angle=90,  margin(0,0,0,0))) 
nn_std_colr <- 'grey70'

get_metric_nn <- function(metric, m_rib) {
  m <- bind_rows(lapply(split(avg_v, avg_v$label), function(i){
    a <- order(i[, metric])
    i$order[a] <- seq(1,num_boots,1)
    i}))
  
  m2 <- left_join(avg_t, m[, c('label','boot','order')], by=c('label','boot'))
  m <- bind_rows(m2, m)
  m <- left_join(m, m_rib)
  m$label <- factor(m$label, levels=c('enhancing', 'nCET', 'necrosis', 'edema',   'expansive', 'focal'))
  m
}

# ROC PLOT
roc_rib <- nn %>% group_by(partition, label, boot) %>% 
  dplyr::summarize(sd=sd(roc)) %>% 
  dplyr::select(partition, label, boot, sd)
roc <- get_metric_nn(metric='roc',m_rib=roc_rib)

nn_roc <- ggplot(roc, aes(y=roc, x=order)) + 
  geom_ribbon(aes(ymin=roc-sd, ymax=roc+sd), fill=nn_std_colr, alpha=0.4) + geom_line(aes(color=label)) + 
  geom_hline(yintercept=0.5, color='black', linetype='dashed') +
  nn_theme + guides(colour=F) +
  # scale_y_continuous(breaks=seq(0.3, 1.0, .1), limits=c(0.3,1.2)) +
  ylab('mean AUC') + xlab('bootstraps') +
  facet_grid(label~partition); print(nn_roc)

# PR PLOT
pr_rib <- nn %>% group_by(partition, label, boot) %>% dplyr::summarize(sd=sd(pr)) %>% dplyr::select(partition, label, boot, sd)
pr <- get_metric_nn(metric='pr',m_rib=pr_rib)

nn_pr <- ggplot(pr, aes(y=pr, x=order)) + 
  geom_ribbon(aes(ymin=pr-sd, ymax=pr+sd), fill=nn_std_colr, alpha=0.4) + geom_line(aes(color=label)) + 
  nn_theme + #guides(colour=F) +
  scale_y_continuous(breaks=seq(0.0, 1.0, .2)) +
  ylab('mean average precision') + xlab('bootstraps') +
  facet_grid(label~partition); print(nn_pr)

# LOSS PLOT
loss_rib <- nn %>% group_by(partition, label, boot) %>% dplyr::summarize(sd=sd(loss)) %>% dplyr::select(partition, label, boot, sd)
loss <- get_metric_nn(metric='loss',m_rib=loss_rib)

nn_loss <- ggplot(loss, aes(y=loss, x=order)) + 
  geom_ribbon(aes(ymin=loss-sd, ymax=loss+sd), fill=nn_std_colr, alpha=0.4) + geom_line(aes(color=label)) + 
  nn_theme + guides(colour=F) +
  # scale_y_continuous(breaks=seq(0.0, 1.0, .2)) +
  ylab('mean loss') + xlab('bootstraps') +
  facet_grid(label~partition); print(nn_loss)

# EPOCH PLOT
ep_rib <- nn %>% group_by(partition, label, boot) %>% dplyr::summarize(sd=sd(epoch)) %>% dplyr::select(partition, label, boot, sd)
ep <- get_metric_nn(metric='epoch', m_rib=ep_rib)

nn_ep <- ggplot(ep[ep$partition=='training',], aes(y=epoch, x=order)) + 
  geom_ribbon(aes(ymin=epoch-sd, ymax=epoch+sd), fill=nn_std_colr, alpha=0.4) + geom_line(aes(color=label)) + 
  nn_theme + guides(colour=F) +
  # scale_y_continuous(breaks=seq(0.0, 1.0, .2)) +
  ylab('mean epoch') + xlab('bootstraps') +
  facet_grid(label~partition); print(nn_ep)

# ACC PLOT
acc_rib <- nn %>% group_by(partition, label, boot) %>% dplyr::summarize(sd=sd(acc)) %>% dplyr::select(partition, label, boot, sd)
acc <- get_metric_nn(metric='acc', m_rib=acc_rib)

nn_acc <- ggplot(acc, aes(y=acc, x=order)) + 
  geom_ribbon(aes(ymin=acc-sd, ymax=acc+sd), fill=nn_std_colr, alpha=0.4 ) + geom_line(aes(color=label)) + 
  nn_theme + guides(fill = guide_legend(title=' a'),
                    color= guide_legend(title='label')) + #  guides(colour=F) +
  # scale_y_continuous(breaks=seq(0.0, 1.0, .2)) +
  ylab('mean accuracy') + xlab('bootstraps') +
  facet_grid(label~partition); print(nn_acc)


plot_grid(nn_ep, nn_loss, nn_roc, nn_pr,
          labels=c('a', 'b','c','d'),
          nrow=1,
          rel_heights=c(.05,1),
          rel_widths=c(.65,1,1,1.2))
dev.print(png, file.path(fdd, paste0('boot_nn_metrics_200','.png')), res=300, height=8, width=13, units="in")


