#' Association tests
#' 
library(dplyr)
library(tidyverse) # combn
library(broom) # map2
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(cowplot) # arrange graphs
library(plyr)
library(foreach)
library(doParallel)

data_dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular/nn_retraining_data'
fdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_figs/association'

# load label info -----
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

# clinical info  --------------
#' 
#'
cdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'
clin <- 'nationwidechildrens.org_clinical_patient_gbm.txt'
clin <- read.table(file.path(cdd, clin), sep='\t', skip=1, header=T, stringsAsFactors=F) # get clinical
clin <- clin[-1,]
clin$X <- clin$bcr_patient_barcode
clin$age_at_initial_pathologic_diagnosis <- as.numeric(clin$age_at_initial_pathologic_diagnosis)
mean_age <- mean(clin$age_at_initial_pathologic_diagnosis)
ma <- as.character(round(mean_age, digit=0))
clin$'diagnosis age' <- unlist(lapply(clin$age_at_initial_pathologic_diagnosis, function(i) ifelse(i < mean_age, 
                                                                                                   paste('<', ma), 
                                                                                                   paste('>=', ma))))
clin$karnofsky_performance_score <- as.numeric(clin$karnofsky_performance_score)
clin$eastern_cancer_oncology_group <- as.numeric(clin$eastern_cancer_oncology_group)


# association data ---------
# test one
table(labels$f6, labels$f10)
fisher.test(labels$f6, labels$f10)

## vasari with vasari
a <- data.frame(t(combn(names(labels),2)), stringsAsFactors = F) %>%
  mutate(d = map2(X1, X2, ~tidy(fisher.test(labels[,.x], labels[,.y])))) %>%
  unnest()
a$p <- p.adjust(a$p.value, method='bonferroni')
a$p <- round(a$p, digits=3)

a$X1 <- mapvalues(a$X1, 
                  from=c('f5','f6','f7','f9','f10', 'f14'), 
                  to=c('enhancing', 'nCET','necrosis','focal','infiltrative','edema'))

a$X2 <- mapvalues(a$X2, 
                  from=c('f6','f7','f9','f10', 'f14', 'subtype'), 
                  to=c('nCET','necrosis','focal','infiltrative','edema', 'subtype'))

a$X1 <- factor(a$X1, levels=c('enhancing', 'nCET','necrosis','focal','infiltrative','edema'))
a$X2 <- factor(a$X2, levels=rev(c('nCET','necrosis','focal','infiltrative','edema', 'subtype')))
a$sig <- ifelse(a$p < 0.05, 'yes','no')

## vasari with clinical
these <- c('race','gender','diagnosis age' , 'karnofsky_performance_score')
pairs <- as.data.frame(expand.grid(names(labels), these, stringsAsFactors=F))
colnames(pairs) <- c('vasari','clinical')

labs <- labels
labs$X <- rownames(labs)
clin_vasari <- left_join(labs, clin[, c('X',these)])
clin_vasari[clin_vasari=='[Not Available]'] <- NA


table(clin_vasari$f7, clin_vasari$'race')
fisher.test(clin_vasari$f7, clin_vasari$'race')

table(clin_vasari$f9, clin_vasari$'karnofsky_performance_score')
fisher.test(clin_vasari$f9, clin_vasari$'karnofsky_performance_score')

table(clin_vasari$f14, clin_vasari$'tissue_source_site')
fisher.test(clin_vasari$f14, clin_vasari$'tissue_source_site')

table(clin_vasari$subtype, clin_vasari$'tissue_source_site')
fisher.test(clin_vasari$subtype, clin_vasari$'tissue_source_site', workspace = 2e8)

c <- pairs %>%
  mutate(d = map2(vasari, clinical, ~tidy(fisher.test(clin_vasari[,.x], clin_vasari[,.y], workspace = 2e8)))) %>%
  unnest()

c$p <- p.adjust(c$p.value, method='bonferroni')
c$p <- round(c$p, digits=3)

c$vasari <- mapvalues(c$vasari, 
                  from=c('f5','f6','f7','f9','f10', 'f14'), 
                  to=c('enhancing', 'nCET','necrosis','focal','infiltrative','edema'))

c$clinical <- mapvalues(c$clinical, 
                      from=c('karnofsky_performance_score','tissue_source_site', 'diagnosis age'), 
                      to=c('Karnofsky\nperformance\nscore', 'tissue\nsource\nsite', 'diagnosis\nage'))

c$sig <- ifelse(c$p < 0.05, 'yes','no')


# plot ---------

fs <- 8
th <- theme_minimal() + 
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        legend.position='bottom',
        axis.text=element_text(size=fs),
        legend.text=element_text(size=fs, color="grey30"),
        legend.title=element_text(size=fs),
        legend.key.size = unit(0.4, "cm"),
        legend.margin=margin(l=-.2, r=0, unit='cm'))

# vasari-vasari ----
# p-values
pv <- ggplot(a, aes(x=X2, y=X1, fill=sig)) + 
  geom_tile(color='grey60') +
  geom_text(aes(label = round(p, 3)), size=3, color='grey30') +
  scale_fill_manual(values = wes_palette("Royal1")[3:4]) +
  labs(fill='adjusted p-value < 0.05') +
  th
  
# odds ratio, for binary only
ev <- ggplot(a, aes(x=X2, y=X1, fill=estimate)) + 
  geom_tile(color='grey60') +
  geom_text(aes(label = round(estimate, 3)), size=3, color='white') +
  scale_fill_gradientn(colors=brewer.pal(n=9, name='RdPu'), breaks=c(0,1,5), limits=c(0,5))+
  th

a_plot <- plot_grid(pv, ev, nrow=1, rel_widths=c(1,1), labels=c('a'), label_size=14)

ff <- file.path(fdd, 'label_association_fisher.png');print(ff)
dev.print(png, ff, res=300, height=3, width=7.5, units="in")


# vasari-clinical  ----
# p-values
pv <- ggplot(c, aes(x=clinical, y=vasari, fill=sig)) + 
  geom_tile(color='grey60') +
  geom_text(aes(label = round(p, 3)), size=3, color='grey30') +
  scale_fill_manual(values = wes_palette("Royal1")[3:4]) +
  labs(fill='adjusted p-value < 0.05') +
  th

# odds ratio, for binary only
ev <- ggplot(c, aes(x=clinical, y=vasari, fill=estimate)) + 
  geom_tile(color='grey60') +
  geom_text(aes(label = round(estimate, 3)), size=3, color='white') +
  scale_fill_gradientn(colors=brewer.pal(n=9, name='RdPu'), breaks=c(0,1,2), limits=c(0,2))+
  th

c_plot <- plot_grid(pv, ev, nrow=1, rel_widths=c(1,1), labels=c('b'), label_size=14)

ff <- file.path(fdd, 'label_clinical_fisher.png');print(ff)
dev.print(png, ff, res=300, height=3, width=7.5, units="in")

# combined ----
plot_grid(a_plot, c_plot, nrow=2, rel_widths=c(1,1))

ff <- file.path(fdd, 'label_combined_fisher.png');print(ff)
dev.print(png, ff, res=300, height=6, width=7.5, units="in")


