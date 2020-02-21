library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(wesanderson)
library(rcartocolor)
library(scales)
library(cowplot) # arrange graphs
library(fgsea)
library(foreach)
library(doParallel)

# dirs
dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular'
fdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_figs/gsea'
data_dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular/nn_retraining_data'
gs_dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'

clean_gene_ids <- function(ids) {
  ids <- gsub("-", ".", ids)
  ids <- gsub("@", ".", ids)
}

rename_classes <- function(c) {
  k <- c
  if (c=='f5') k <- 'enhancing'
  if (c=='f6') k <- 'nCET'
  if (c=='f7') k <- 'necrosis'
  if (c=='f9') k <- 'focal'
  if (c=='f10') k <- 'infiltrative'
  if (c=='f14') k <- 'edema'
  return(k)
}

rename_gene_sets <- function(names, type) {
  n <- names
  if (!(type %in% c('verhaak','single'))) {
    n <- tolower(n)
    n <- gsub('hallmark_','', n)
    n <- gsub('reactome_','', n)
    n <- gsub('go_','', n)
    n <- gsub('_',' ', n)
    n <- gsub('^barres', 'Zhang', n)
    n <- gsub('^darmanisbarres', 'Darmanis', n)
    n <- gsub('^patel', 'Patel', n)
    n <- gsub('gbmcoreastrocytes', 'gbm core astrocytes', n)
    n <- gsub('matureastrocytes', 'mature astrocytes', n)
    n <- gsub('microgliamacrophages', 'microglia macrophages', n)
  }
  return(n)
}

# get literature related findings -------
## radiogenomics

## function to read saved query results
get_search <- function(fname, tag) {
  print(fname)
  gs <- gmtPathways(gmt.file=file.path(dd,'model_ready','msigdb_v6.2_GMTs','search_results', paste0(fname,'.gmt')))
  cat('num gs: ', length(gs))
  gs <- data.frame(X=na.omit(names(gs)), query=fname)
}

diehn_enhan <- c('hypoxia', 'ECM','angiogenesis')
diehn_infil <- c('immune')
gev_necro <- c('GAP43', 'IL4', 'cell_membrane') # gevaert

# jamshidi
jam_enhan_paths <- c('BIOCARTA_AGR_PATHWAY','BIOCARTA_AT1R_PATHWAY', 'BIOCARTA_BIOPEPTIDES_PATHWAY','BIOCARTA_CARDIACEGF_PATHWAY','BIOCARTA_CREB_PATHWAY',
                     'BIOCARTA_EGF_PATHWAY','BIOCARTA_ERK5_PATHWAY','BIOCARTA_IGF1R_PATHWAY','BIOCARTA_IL7_PATHWAY','BIOCARTA_INTRINSIC_PATHWAY',
                     'BIOCARTA_MAPK_PATHWAY', 'BIOCARTA_MET_PATHWAY','BIOCARTA_NFAT_PATHWAY','BIOCARTA_PDGF_PATHWAY', 'BIOCARTA_TFF_PATHWAY', 
                     'BIOCARTA_TOLL_PATHWAY')
jam_enhan_genes <- c('C1orf172', 'CAMSAP2', 'KCNK3', 'LTBP1')
jam_necro_paths <- c('BIOCARTA_CK1_PATHWAY', 'BIOCARTA_CYTOKINE_PATHWAY','BIOCARTA_PML_PATHWAY','BIOCARTA_TID_PATHWAY')
jam_necro_genes <- c('ITGA5', 'RUNX3')

gut_enhan <- c('EGFR') # gutman
gut_necro <- c('CDKN2A')

colen_invas <- c('MYC', 'NFKBIA')
zinn_invas_up <- c( 'POSTN', 'CXCL12', 'COL1A1', 'COL6A3', 'GRB10')
zinn_invas_dn <- c( 'AQP4', 'KIF1A', 'MPPED2', 'CTNNA2', 'PKP4', 'KIF5C')

## genomics
parsons <- c('TP53', 'PI3K','RB1','EGFR','PTEN', 'NF1','CDK4','IDH1','PIK3CA', 'PIK3R1', 'MGMT_promoter') # added mgmt myself

# patient info  --------------
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

# load patient label info
lab_fns <- c('f5','f6','f7','f9','f10','f14','subtype')
parallel::mcaffinity(1:7)
cl <- makeCluster(7)
registerDoParallel(cl)
labels <-  foreach (i=1:length(lab_fns), .export=c()) %dopar% {
  read.csv(file.path(data_dd, paste0(lab_fns[i],'.csv')), stringsAsFactors=F)[, 1:2]
}
stopCluster(cl)

a <- Reduce(function(x,y) full_join(x, y, by='X'), labels)
colnames(a) <- lapply(colnames(a), function(i) strsplit(i, '\\.')[[1]][1] )
these <- c('X','race','gender','ethnicity','diagnosis age' , 'tissue_source_site')
clin <- full_join(a, clin[, these], by='X')
rm(a)
clin$X <- gsub('-', '.' ,clin$X)


#  gene set info  --------------
#' 
#'
gs_dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'

verhaak <- read.table(file.path(gs_dd, 'TCGA_unified_CORE_ClaNC840.txt'), sep='\t', header=T, stringsAsFactors=F)
verhaak <- verhaak[-1,1:2]
colnames(verhaak) <- c('go_genes', 'subtype')
verhaak$go_genes <- clean_gene_ids(as.character(verhaak$go_genes))
verhaak$X <- verhaak$go_genes
t <- lapply(c('PN','NL','CL','MES'), function(i) verhaak$go_genes[verhaak$subtype==i])
t[['all']] <- verhaak$go_genes
verhaak <- t
names(verhaak) <- c('PN','NL','CL','MES', 'all')

puchalski <- read.table(file.path(gs_dd, 'gene_sets_Puchalski', 'aaf2666_Table-S15.csv'), sep=',', header=T, stringsAsFactors=F)
puchalski <- puchalski[-1,]
puchalski <- lapply(puchalski, clean_gene_ids)
puchalski <- as.data.frame(puchalski, stringsAsFactors=F)
# BI gene sets

gs_files <- list('h.all.v6.2.symbols.gmt','c1.all.v6.2.symbols.gmt',
                 'c2.cp.reactome.v6.2.symbols.gmt', 'c2.cp.biocarta.v6.2.symbols.gmt','c2.cp.kegg.v6.2.symbols.gmt',
                 'c2.cp.v6.2.symbols.gmt', 'c2.cgp.v6.2.symbols.gmt',
                 'c3.all.v6.2.symbols.gmt',
                 'c4.all.v6.2.symbols.gmt', 'c5.all.v6.2.symbols.gmt', 'c6.all.v6.2.symbols.gmt', 'c7.all.v6.2.symbols.gmt')
names(gs_files) <- c('hallmark','chromosome',
                     'reactome', 'biocarta','kegg',
                     'canonical', 'chem',
                     'motif', 'computational', 'GO','onco','immuno')

# get gene ids 
gids <- clean_gene_ids(read.csv(file.path(dd,'nn_gene_masking', 'f5_single_scores.csv'), stringsAsFactors=F)$X)

gs_names <- c(names(gs_files), 'verhaak','puchalski')
gs_info <- lapply(gs_names, function(n) {
  print(n)
  if (n=='verhaak') {
    subs <- names(verhaak)
    g <- lapply(subs, function(i) {
      gs <- verhaak[[i]]
      gs[!(gs=='')]
    })
    names(g) <- subs
  } else if (n=='puchalski') {
    g <- lapply(colnames(puchalski), function(i) {
      gs <- puchalski[,i]
      gs[!(gs=='')]
    })
    names(g) <- colnames(puchalski)
  } else {
    g <- gmtPathways(gmt.file= file.path(gs_dd,'msigdb_v6.2_GMTs', gs_files[[n]]) )
  }
  
  size <- unlist(lapply(g, function(i) length(i))) # number of genes in each gene set 
  cover_n <- unlist(lapply(g, function(i) length( i[i %in% gids] ) )) # number of genes in gene set found in microarray data
  cover_p <- round(100*cover_n/size, digits=1)
  
  genes <- unique(unlist(g))  
  gs_size <- length(genes) # total genes in overall gene set
  na_genes <-  genes[!(genes %in% gids)]   # genes not found in micrroary
  na_size <- length(na_genes)
  
  inf <- data.frame(cover=cover_p, 'size'=size)
  inf$'gene set' <- n
  inf$X <- names(g)
  if (n=='verhaak') {
    rand_info <- data.frame(cover=100, size=200, X='random') # attach the 200 random, non-verhaak genes info
    rand_info$'gene set' <- 'verhaak'
    inf <- rbind(inf, rand_info)
  }
  
  stats <- list('total'=gs_size,'na_genes'=na_genes, 'na_size'=na_size )
  list(info=inf, stats=stats)
})
names(gs_info) <- gs_names

gs_stats <- lapply(gs_info, '[[', 2)
gs_info <- lapply(gs_info, '[[', 1)
gs_info <- bind_rows(gs_info)

# # remove duplicates
# gs_info <- split(gs_info, gs_info$`gene set`)
# gs_info$canonical <- gs_info$canonical[ !(gs_info$canonical$X %in% c(gs_info$biocarta$X,gs_info$kegg$X, gs_info$reactome$X)), ]
# gs_info <- do.call('rbind', gs_info)
# rownames(gs_info) <- gs_info$X


# heatmap colors  ----------
subtype_color <- c('white','grey','red','yellow','blue','green')
names(subtype_color) <- c('na','unnamed','NL','PN','MES','CL')

gene_set_color <- c(RColorBrewer::brewer.pal(name='Paired',n=12), RColorBrewer::brewer.pal(name='Set3',n=4))
names(gene_set_color) <- c('hallmark','chromosome','motif','computational','GO',
                           'onco','immuno', 'verhaak', 'puchalski', 'single',
                           'reactome', 'random',
                           'kegg','biocarta','canonical','CGP')
unique(gs_info$`gene set`) %in% names(gene_set_color)



# theme ----
t <- 10; t2 <- t + 0
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
# saliency gsea ----------
these_gs <- c('verhaak','puchalski', 'hallmark','reactome', 'chromosome', 'GO', 'motif','onco','immuno', 'computational')
# these_gs <- c('verhaak','puchalski', 'hallmark','reactome', 'chromosome', 'GO', 'motif')
these_gs <- gs_names

get_heatmap_sal <- function(s, cut_col=4, fs=8, all_annotations=T) {
  r1 <- as.data.frame(s[, c('pathway', 'padj','name', 'class')])
  r1 <- spread(r1, pathway, padj)
  rownames(r1) <- r1$name
  
  r2 <- as.data.frame(s[, c('pathway', 'NES','name', 'class')])
  r2 <- spread(r2, pathway, NES)
  rownames(r2) <- r2$name
  
  all.equal(rownames(r1), rownames(r2))
  
  # save heat annotations
  
  if (r1$class[1]=='') {
    demo <- c(demo, 'subtype', label)
    info <- r1[, c('name'), drop=F] 
    colnames(info) <- c('X')
  } else {
    info <- r1[, c('name','class'), drop=F] 
    colnames(info) <- c('X', label)
  }
  
  
  if (all_annotations) {
    demo <- c('X','diagnosis age','gender', 'ethnicity','race')
    info <- left_join(info, clin[, demo], by='X')
    # all.equal(info[, paste0(label, '.x')], info[, paste0(label, '.y')]) # check matching worked
    # colnames(info)[colnames(info) %in% paste0(label, '.x')] <- label # remove repeated label info
    # info <- info[,colnames(info) != paste0(label, '.y') ]
    
    rownames(info) <- info$X
    info <- info[ , colnames(info) != 'X', drop=F]
    info[is.na(info)] <- '[Not Available]'
    info[info=='[Not Available]'] <- 'na'
  } else {
    rownames(info) <- info$X
    info <- info[ , colnames(info) != 'X', drop=F]
  }


  
  # heatmap annotation colors
  ann_colors <- list('subtype'=c('na'='grey','Neural'='red','Proneural'='yellow','Mesenchymal'='blue','Classical'='green'), 
                     'diagnosis age'=c('na'='grey','< 58'='#fdbf6f','>= 58'='#ff7f00'),
                     'gender'=c('na'='grey', 'FEMALE'='#fb9a99','MALE'='#a6cee3'),
                     'ethnicity'=c('na'='grey','HISPANIC OR LATINO'='#b15928','NOT HISPANIC OR LATINO'='#ffff99'),
                     'race'=c('na'='grey','WHITE'='#fb9a99', 'BLACK OR AFRICAN AMERICAN'='#e31a1c', 'ASIAN'='#cab2d6'),
                     'f5'=c('na'='grey',' More than 1/3'='#238b45', 'Less than 1/3'='#99d8c9'),
                     'f6'=c('na'='grey',' More than 1/3'='#238b45', 'Less than 1/3'='#99d8c9'),
                     'f7'=c('na'='grey',' More than 1/3'='#238b45', 'Less than 1/3'='#99d8c9'),
                     'f14'=c('na'='grey',' More than 1/3'='#238b45', 'Less than 1/3'='#99d8c9'),
                     'f9'=c('na'='grey','Focal'='#a6bddb', 'non-focal'='#014636'),
                     'f10'=c('na'='grey','exp'='#fde0dd', 'mix'='#f768a1', 
                             'inf'='#7a0177')
  )
  
  # heat mat pre
  r1 <- r1[, !(colnames(r1) %in% c('name','class')) ]  # remove metadata
  r2 <- r2[, !(colnames(r2) %in% c('name','class')) ]
  r1 <- r1[, apply(r1, 2, function(r) any(r < 0.05))]  # keep any cols have sig results
  r2 <- r2[, colnames(r1)]
  r1 <- apply(r1, 2, function(i) {  # binarize significance
    p <- i
    p[i >= 0.05] <- 0
    p[i < 0.05] <-  1
    p
  })
  r1 <- t(r1) # transpose
  r2 <- t(r2)
  
  breaks <- seq(0,2.5,0.1)
  col <- colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(length(breaks))
  p1 <- pheatmap(mat=as.matrix(r2), show_colnames=F, main='normalized enrichment score',
                 cutree_cols=cut_col, treeheight_row=0, breaks=breaks,
                 color=col,
                 legend=T, fontsize=fs, fontsize_row=fs+2, fontsize_col=fs+2,
                 legend_breaks=c(seq(0,3, .5)), legend_labels=c(seq(0,3, .5)),
                 annotation_col=info, annotation_legend=T, annotation_colors=ann_colors)
  
  r11 <- r1[p1$tree_row$order, p1$tree_col$order]
  rownames(r11) <- rownames(r1)[p1$tree_row$order]
  colnames(r11) <- colnames(r1)[p1$tree_col$order]
  
  cls <- cutree(p1$tree_col, cut_col)
  cls <- cls[p1$tree_col$order]
  gaps <- unlist(lapply(unique(cls), function(i) length(cls[cls==i])))
  gaps <- unlist(lapply(seq(1,cut_col, 1), function(i) sum(gaps[1:i])))
                 
  col2 <- colorRampPalette(brewer.pal(n = 5, name = "RdBu"))(2)
  p2 <- pheatmap(mat=as.matrix(r11), show_colnames=F,  main='adjusted p-value',
                 gaps_col=gaps, cluster_rows=F, cluster_cols=F,
                 color=rev(col2),
                 legend=T, fontsize=fs, fontsize_row=fs+2, fontsize_col=fs+2,
                 legend_breaks=c(0,1), legend_labels=c('p >= 0.05','p < 0.05'),
                 annotation_col=info[, label, drop=F], annotation_legend=F, annotation_colors=ann_colors)
  return(list(p1,p2))
}

# plot gsea heatmap -----------------------
these_gs <- c('verhaak','hallmark','puchalski','reactome')
fs_gs <- c(10, 7, 8, 6)
names(fs_gs) <- these_gs

label <- 'subtype' 
subtype_sal <- lapply(these_gs, function(gs_name){
  f <- file.path(dd,'gsea', paste('saliency', label, gs_name, '.rds', sep='_')); print(f)
  s <- readRDS(f) 
  get_heatmap_sal(s=s$sal, fs=fs_gs[gs_name])
})
names(subtype_sal) <- these_gs

dev.new();
gs_name <- 'verhaak'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
a <- plot_grid(subtype_sal[[gs_name]][[2]][[4]],NULL, nrow=1, rel_widths=c(1, 0.225))
plot_grid(subtype_sal[[gs_name]][[1]][[4]], a, nrow=2, rel_heights=c(1, .7))
dev.print(png, ff, res=300, height=6, width=18, units="in")

gs_name <- 'puchalski'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
a <- plot_grid(NULL,  subtype_sal[[gs_name]][[2]][[4]], NULL, nrow=1, rel_widths=c(0.045, 1, 0.225))
plot_grid(subtype_sal[[gs_name]][[1]][[4]], a, nrow=2, rel_heights=c(1, .7))
dev.print(png, ff, res=300, height=6, width=18, units="in")
# 
# gs_name <- 'hallmark'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
# a <- plot_grid(NULL, subtype_sal[[gs_name]][[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.225))
# plot_grid(subtype_sal[[gs_name]][[1]][[4]], a, nrow=2, rel_heights=c(1, .85))
# dev.print(png, ff, res=300, height=18, width=20, units="in")
# 
# gs_name <- 'reactome'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
# a <- plot_grid(NULL, subtype_sal[[gs_name]][[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.125))
# plot_grid(psubtype[[gs_name]][[4]], a, nrow=2, rel_heights=c(1, .85))
# dev.print(png, ff, res=300, height=25, width=25, units="in")
# 
# gs_name <- 'GO'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
# a <- plot_grid(NULL, subtype_sal[[gs_name]][[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.125))
# plot_grid(subtype_sal[[gs_name]][[1]][[4]], a, nrow=2, rel_heights=c(1, .85))
# dev.print(png, ff, res=300, height=30, width=25, units="in")
#
label <- 'f5'
f5 <- lapply(these_gs, function(gs_name){
  f <- file.path(dd,'gsea', paste('saliency', label, gs_name, '.rds', sep='_')); print(f)
  s <- readRDS(f) 
  get_heatmap_sal(s=s$sal, title=gs_name)
})
names(f5) <- these_gs

gs_name <- 'verhaak'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
p <- f5[[gs_name]]
a <- plot_grid(NULL, p[[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.225))
plot_grid(p[[1]][[4]], a, nrow=2, rel_heights=c(1, .7))
dev.print(png, ff, res=300, height=7, width=18, units="in")

gs_name <- 'puchalski'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
p <- f5[[gs_name]]
a <- plot_grid(NULL, p[[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.225))
plot_grid(p[[1]][[4]], a, nrow=2, rel_heights=c(1, .7))
dev.print(png, ff, res=300, height=7, width=18, units="in")

# gs_name <- 'hallmark'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
# p <- f5[[gs_name]]
# a <- plot_grid(NULL, p[[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.225))
# dev.new(); plot_grid(p[[1]][[4]], a, nrow=2, rel_heights=c(1, .85))
# dev.print(png, ff, res=300, height=12, width=18, units="in")
# 
# gs_name <- 'reactome'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
# p <- f5[[gs_name]]
# a <- plot_grid(NULL, p[[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.125))
# plot_grid(p[[1]][[4]], a, nrow=2, rel_heights=c(1, .85))
# dev.print(png, ff, res=300, height=25, width=25, units="in")
# 
# gs_name <- 'GO'; ff <- file.path(fdd, paste('saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
# p <- f5[[gs_name]]
# a <- plot_grid(NULL, p[[2]][[4]],NULL, nrow=1, rel_widths=c(0.045, 1, 0.125))
# dev.new();plot_grid(p[[1]][[4]], a, nrow=2, rel_heights=c(1, .85))
# dev.print(png, ff, res=300, height=30, width=25, units="in")




# get patient counts of enriched gs for each label ----------

get_counts <- function(gs_name, count_thres=5) {
  vasari_labs <- c('f5','f6','f7','f9','f10','f14')
  # cols <- carto_pal(9, 'Bold')[c(2,3,1,6,5,7,4,8,9)]
  cols <- carto_pal(12, 'Prism')[c(2,3,5,6,8,9,10,12)]
  
  counts <- lapply(vasari_labs, function(lab) {
    f <- file.path(dd,'gsea', paste('saliency', lab, gs_name, '.rds', sep='_')); print(f)
    s <- readRDS(f)
    
    c <- s$sal[s$sal$padj < 0.05,] 
    if (nrow(c) > 0) {
      c <- c %>% dplyr::group_by(pathway) %>% 
        dplyr::summarize(count=length(name)) %>%
        dplyr::select(pathway, count)
      c$label <- lab
    } else {
      print('no significant associations')
      return(c)
    }
    c
  })
  counts <- bind_rows(counts)
  if (nrow(counts) > 0)
    counts$name <- gs_name
    
    counts$pathway <- rename_gene_sets(counts$pathway, type=gs_name)
    counts$pathway <- gsub(' ', '\n', counts$pathway)
    
    counts$label <- unlist(lapply(counts$label, rename_classes))
    counts$label <- factor(counts$label, levels=c('enhancing','nCET', 'necrosis','edema','focal','infiltrative')) # order to plot consistent colors
    these_gs <- unique(counts$pathway[counts$count >= count_thres])
    counts <- counts[counts$pathway %in% these_gs, ]
    
    p <- ggplot(counts, aes(x=pathway, y=count, fill=label)) +
      geom_bar(stat='identity', position=position_dodge2(width=0.9, preserve='single', padding=0), color = "black", lwd=0.3) +
      # geom_text(aes(label=count, x=pathway),  position=position_dodge2(width=0.9, preserve='single'), vjust=-.1, color='gray30') +
      theme_bw() + th + 
      # coord_flip() +
      scale_fill_manual(values=cols) +
      labs(x='gene set', title=gs_name); print(p)
    
    ff <- file.path(fdd, paste('saliency', 'count', 'thres', count_thres, gs_name, 'barplot.png', sep='_')); print(ff)
    dev.print(png, ff, res=300, height=3, width=9, units="in")
  return(counts)
}


names(gs_files) # available

these_gs <- c('verhaak','puchalski', names(gs_files))
thresholds <- c(rep(10, length(these_gs)))
names(thresholds) <- c(these_gs)
thresholds['chem'] <- 20
thresholds['computational'] <- 20
# thresholds['chromosome'] <- 5
# thresholds['hallmark'] <- 5
# thresholds['immuno'] <- 5
# thresholds['onco'] <- 5
# thresholds['reactome'] <- 10

z <- lapply(names(thresholds) , function(i) {
  print(i)
  get_counts(gs_name=i, count_thres=thresholds[i])
})
plotz <- bind_rows(z)


a <- get_counts(gs_name='puchalski', count_thres=1)
b <- get_counts(gs_name='hallmark', count_thres=1)

get_counts(gs_name='chromosome', count_thres=10)
get_counts(gs_name='canonical', count_thres=5)


## combine multiple count plots
unique(plotz$name)
plotz$name <- gsub('verhaak','subtype', plotz$name)
plotz$name <- gsub('puchalski','cell types or phenotypes', plotz$name)
plotz$pathway <-  gsub('fetalneuronsreplicating','fetal\nneurons\nreplicating', plotz$pathway)
plotz$pathway <-  gsub('anticellcycle','anticell\ncycle', plotz$pathway)
plotz$pathway <-  gsub('oligodendrocytes','oligo-\ndendrocytes', plotz$pathway)
plotz$pathway <-  gsub('checkpoint','check-\npoint', plotz$pathway)
plotz$pathway <-  gsub('trafficking\nof\nampa\nreceptors','trafficking of ampa receptors', plotz$pathway)

plotz$name <- factor(plotz$name, levels=c('subtype','cell types or phenotypes','hallmark','chromosome','reactome','onco','immuno', 'chem', 'computational'))
cols <- carto_pal(12, 'Prism')[c(2,5,6,8,3,11)]
cols <- carto_pal(12, 'Prism')[c(5,6,8,3,11)]



# 'reactome','onco', 'chem', 'computational' do in supp
ggplot(plotz[plotz$name %in% c('subtype','cell types or phenotypes','hallmark'), ], 
       aes(x=pathway, y=count, fill=label)) +
  geom_bar(stat='identity', position=position_dodge2(width=0.9, preserve='single', padding=0), color = "black", lwd=0.3) +
  geom_text(aes(label=count, x=pathway),  position=position_dodge2(width=0.9, preserve='single'), vjust=-.1, color='gray30') +
  theme_bw() + th + theme(legend.position='right',
                          axis.title.x=element_text(size=t)) +
  scale_fill_manual(values=cols) +
  labs(x='gene sets',
       y='#\nenriched\npatients',
       fill='radiogenomic\nmodel') + 
  # guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  # coord_flip() +
  facet_wrap(~name, scale='free_x', ncol=1)

ff <- file.path(fdd, paste('saliency', 'percent', 'vasari', 'barplot_byclass_text.png', sep='_')); print(ff)
ff <- file.path(fdd, paste('saliency', 'percent', 'vasari', 'barplot_byclass.png', sep='_')); print(ff)
dev.print(png, ff, res=300, height=7, width=9, units="in")

# subtype
lab <- 'subtype'
gs_name <- 'verhaak'
f <- file.path(dd,'gsea', paste('saliency', lab, gs_name, '.rds', sep='_')); print(f)
s <- readRDS(f)
c <- s$sal[s$sal$padj < 0.05,] %>% group_by(pathway, class) %>% 
  summarize(count=length(name)) %>%
  dplyr::select(pathway, count, percent, class)

cols <- carto_pal(9, 'Bold')[c(2,3,5,4)]
cplot <-  ggplot(c, aes(x=class, y=count, fill=pathway)) +
    geom_bar(stat='identity', position=position_dodge2(width=0.9, preserve='single', padding=0), color = "black", lwd=0.3) + 
  # geom_text(aes(label=count, x=class),  position=position_dodge2(width=0.9, preserve='single'), vjust=-.1, color='gray30') +
    # coord_flip() +
    scale_fill_manual(values=cols) +
    theme_bw() + th + theme(axis.title.x=element_text(size=t)) +
    labs(x="patients' true subtype",
         y='#\nenriched\npatients',
         fill='gene set'); print(cplot)
ff <- file.path(fdd, paste('saliency', 'percent', 'subtype', gs_name, 'barplot_byclass.png', sep='_')); print(ff)
ff <- file.path(fdd, paste('saliency', 'percent', 'subtype', gs_name, 'barplot_byclass_text.png', sep='_')); print(ff)
dev.print(png, ff, res=300, height=1.5, width=6, units="in")

# figure publish subtype model verhaak gene set saliency ----

these_gs <- c('verhaak')
fs_gs <- c(8)
names(fs_gs) <- these_gs

label <- 'subtype' 
subtype_sal <- lapply(these_gs, function(gs_name){
  f <- file.path(dd,'gsea', paste('saliency', label, gs_name, '.rds', sep='_')); print(f)
  s <- readRDS(f) 
  get_heatmap_sal(s=s$sal, fs=fs_gs[gs_name], all_annotations=F)
})
names(subtype_sal) <- these_gs

gs_name <- 'verhaak'; ff <- file.path(fdd, paste('combined_saliency', label, gs_name, 'heat.png', sep='_')); print(ff)
a <- plot_grid(subtype_sal[[gs_name]][[2]][[4]],NULL, nrow=1, rel_widths=c(1, 0.075))
p1 <- plot_grid(subtype_sal[[gs_name]][[1]][[4]], a, nrow=2, rel_heights=c(1, .7)); print(p1)

p2 <- plot_grid(cplot, NULL, nrow=1, rel_widths=c(1,0.2))
plot_grid(p1, p2, nrow=2, rel_heights=c(1, .4), labels=c('a','b'), label_size=16, scale=0.97)
dev.print(png, ff, res=300, height=5, width=8, units="in")
# 

# figure publish vasari model gene set saliency counts ----
ff <- file.path(fdd, paste('combined_saliency', label, gs_name, 'heat.png', sep='_')); print(ff)

dev.print(png, ff, res=300, height=5, width=8, units="in")
# 

# single and corr gsea ----------
get_heatmap <- function(s, fs=7, cheight=10, rheight=10, row_order=NA, col_order=NA, lim=6) {
  s$class <- lapply(s$class, rename_classes)
  r1 <- as.data.frame(s[, c('pathway', 'padj','class', 'name')])
  r1 <- spread(r1, pathway, padj)
  rownames(r1) <- r1$class

  r2 <- as.data.frame(s[, c('pathway', 'NES','class', 'name')])
  r2 <- spread(r2, pathway, NES)
  rownames(r2) <- r2$class

  all.equal(rownames(r1), rownames(r2))
  
  # heat mat pre
  r1 <- r1[, !(colnames(r1) %in% c('name','class')) ]  # remove metadata
  r2 <- r2[, !(colnames(r2) %in% c('name','class')) ]
  r1 <- r1[, apply(r1, 2, function(r) any(r < 0.05))]  # keep any cols have sig results
  r2 <- r2[, colnames(r1)]
  r1 <- apply(r1, 2, function(i) {  # binarize significance
    p <- i
    p[i >= 0.05] <- 0
    p[i < 0.05] <-  1
    p
  })
  r1 <- t(r1) # transpose
  r2 <- t(r2)
  
  if (!(is.na(row_order)) & !(is.na(col_order))) {
    r1 <- r1[row_order, col_order]
    r2 <- r2[row_order, col_order]
    clus_rows <- F
    clus_cols <- F
  } else {
    clus_rows <- T
    clus_cols <- T
  }
  
  # heat info
  info <- data.frame('X'=rownames(r1))
  info[,'gene set'] <- s$name[1]
  info <- left_join(info, gs_info, by=c('X','gene set'))
  rownames(info) <- rownames(r1) # left join removes rownames, add them back
  info$cover_new <- cut(info$cover, breaks=seq(0,100,10))
  maxsize <- max(na.omit(info$size))
  interval <- 100
  if (maxsize <= 100) interval <- 10
  info$size_new <- cut(info$size,   breaks=seq(0,maxsize+interval, interval))
  
  
  # heat colors
  nc <- levels(droplevels(info$cover_new)) 
  gene_cover_color <- rev(colorRampPalette(rev(brewer.pal(n=6, name='Greens')))(length(nc)))
  names(gene_cover_color) <- nc

  ns <- levels(droplevels(info$size_new)) 
  gene_size_color <- colorRampPalette(brewer.pal(n=6, name="Blues" ))(length(ns))
  names(gene_size_color) <- ns
  
  ann_row <- info[, rev(c('gene set', 'size_new','cover_new'))]
  rownames(ann_row) <- rownames(info)
  ann_row$'% coverage' <- as.character(ann_row$cover_new)
  ann_row$'num. genes' <- as.character(ann_row$size_new)
  ann_row <- ann_row[, rev(c('gene set', '% coverage', 'num. genes'))]
  ann_colors <- list('gene set'=gene_set_color[names(gene_set_color) %in% info$'gene set'],
                     '% coverage'=gene_cover_color,
                     'num. genes'=gene_size_color )
  
  rn <- info$'gene set' %in% c('verhaak', 'single', 'puchalski') # rename for plotting
  lr <- rownames(info)
  lr[!(rn)] <- tolower(lr[!(rn)]) 
  lr <- gsub('hallmark_','', lr)
  lr <- gsub('reactome_','', lr)
  lr <- gsub('go_','', lr)
  lr <- gsub('_',' ', lr)
  
  if (info$'gene set'[1]=='puchalski') {
    # match naming in paper
    lr <- gsub('^Barres', 'Zhang', lr)
    lr <- gsub('^DarmanisBarres', 'Darmanis', lr)
  }

  breaks <- seq(-lim,lim,1)
  col <- colorRampPalette(brewer.pal(n = 11, name = "PiYG"))(length(breaks))
  p1 <- pheatmap(mat=as.matrix(r2), main='normalized enrichment score',
                 cluster_rows=clus_rows, cluster_cols=clus_cols,
                 color=col, breaks=breaks,
                 labels_row=lr, 
                 legend=T, fontsize=fs, fontsize_row=fs+2, fontsize_col=fs+2,
                 treeheight_row=rheight, treeheight_col=0,
                 annotation_row=ann_row, annotation_colors=ann_colors)
  

  
  col2 <- colorRampPalette(brewer.pal(n = 5, name = "RdBu"))(2)
  
  if (!(is.na(row_order)) & !(is.na(col_order))) {
    r11 <- r1
    
    lr2 <- lr
  } else {
    r11 <- r1[p1$tree_row$order, p1$tree_col$order]
    rownames(r11) <- rownames(r1)[p1$tree_row$order]
    colnames(r11) <- colnames(r1)[p1$tree_col$order]
    
    lr2 <- lr[p1$tree_row$order]
    
  }
  
  p2 <- pheatmap(mat=as.matrix(r11), main='adjusted p-value',
                 cluster_rows=F, cluster_cols=F,
                 color=rev(col2),
                 labels_row=lr2, 
                 legend=T, fontsize=fs, fontsize_row=fs+2, fontsize_col=fs+2,
                 legend_breaks=c(0,1), legend_labels=c('p >= 0.05','p < 0.05'))
  return(list(p1,p2))
}

these_gs <- c('verhaak','hallmark','puchalski','chromosome','onco', 'motif', 'canonical')
fs_gs <- c(rep(8, length((these_gs))))
names(fs_gs) <- these_gs

label <- 'subtype' # single --------
subtype_single <- lapply(these_gs, function(gn){
  f <- file.path(dd,'gsea', paste('single', label, gn, '.rds', sep='_')); print(f)
  s <- as.data.frame(readRDS(f))
  if (gn=='verhaak') {
    get_heatmap(s=s, fs=fs_gs[gn], row_order=c('NL','PN','CL', 'MES'), col_order=c('Neural','Proneural','Classical','Mesenchymal'))
  } else {
    get_heatmap(s=s, fs=fs_gs[gn])
  }
})
names(subtype_single) <- these_gs

dev.new();

gs_name <- 'verhaak'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_single[[gs_name]][[1]][[4]], subtype_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .7), scale=.97)
dev.print(png, ff, res=300, height=2.5, width=7, units="in")

gs_name <- 'hallmark'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_single[[gs_name]][[1]][[4]], subtype_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), scale=.95)
dev.print(png, ff, res=300, height=6, width=10, units="in")

gs_name <- 'puchalski'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_single[[gs_name]][[1]][[4]], subtype_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), scale=.95)
dev.print(png, ff, res=300, height=4, width=10.5, units="in")

gs_name <- 'reactome'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_single[[gs_name]][[1]][[4]], subtype_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), scale=.95)
dev.print(png, ff, res=300, height=8, width=20, units="in")


label <- 'vasari' # single --------
vasari_single <- lapply(these_gs, function(gn){
  f <- file.path(dd,'gsea', paste('single', label, gn, '.rds', sep='_')); print(f)
  s <- as.data.frame(readRDS(f))
  max_enrich <- max(abs(c(max(s$NES, na.rm=T), min(s$NES, na.rm=T))))
  get_heatmap(s=s, fs=fs_gs[gn], lim=ceiling((max_enrich)))
})
names(vasari_single) <- these_gs

dev.new();

gs_name <- 'verhaak'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_single[[gs_name]][[1]][[4]], vasari_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .7))
dev.print(png, ff, res=300, height=4, width=8, units="in")

gs_name <- 'hallmark'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_single[[gs_name]][[1]][[4]], vasari_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8))
dev.print(png, ff, res=300, height=4, width=10, units="in")

gs_name <- 'puchalski'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_single[[gs_name]][[1]][[4]], vasari_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8))
dev.print(png, ff, res=300, height=3, width=10.5, units="in")

gs_name <- 'chromosome'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_single[[gs_name]][[1]][[4]], vasari_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8))
dev.print(png, ff, res=300, height=5, width=8, units="in")

gs_name <- 'onco'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_single[[gs_name]][[1]][[4]], vasari_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8))
dev.print(png, ff, res=300, height=5, width=11, units="in")

gs_name <- 'motif'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_single[[gs_name]][[1]][[4]], vasari_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8))
dev.print(png, ff, res=300, height=6, width=12, units="in")

gs_name <- 'canonical'; ff <- file.path(fdd, paste('single', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_single[[gs_name]][[1]][[4]], vasari_single[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8))
dev.print(png, ff, res=300, height=5, width=14, units="in")


label <- 'subtype' # correlation -------
subtype_corr <- lapply(these_gs, function(gn){
  f <- file.path(dd,'gsea', paste('corr', label, gn, '.rds', sep='_')); print(f)
  s <- as.data.frame(readRDS(f))
  get_heatmap(s=s, fs=fs_gs[gn])
})
names(subtype_corr) <- these_gs

dev.new()

gs_name <- 'verhaak'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_corr[[gs_name]][[1]][[4]], subtype_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .7), scale=.97)
dev.print(png, ff, res=300, height=2.5, width=7, units="in")

gs_name <- 'hallmark'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_corr[[gs_name]][[1]][[4]], subtype_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), scale=.95)
dev.print(png, ff, res=300, height=6, width=10, units="in")

gs_name <- 'puchalski'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_corr[[gs_name]][[1]][[4]], subtype_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), scale=.95)
dev.print(png, ff, res=300, height=4, width=10.5, units="in")

gs_name <- 'reactome'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(subtype_corr[[gs_name]][[1]][[4]], subtype_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), scale=.95)
dev.print(png, ff, res=300, height=8, width=20, units="in")



label <- 'vasari' # correlation ----------
these_gs <- c('verhaak','hallmark','puchalski','chromosome','onco')
fs_gs <- c(rep(9, length((these_gs))))
names(fs_gs) <- these_gs

vasari_corr <- lapply(these_gs, function(gn){
  f <- file.path(dd,'gsea', paste('corr', label, gn, '.rds', sep='_')); print(f)
  s <- as.data.frame(readRDS(f))
  max_enrich <- max(abs(c(max(s$NES, na.rm=T), min(s$NES, na.rm=T))))
  get_heatmap(s=s, fs=fs_gs[gn], lim=ceiling((max_enrich)))
})
names(vasari_corr) <- these_gs

dev.new();

gs_name <- 'verhaak'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_corr[[gs_name]][[1]][[4]], vasari_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .7), labels=c('a','b'), label_size=16, scale=.97)
dev.print(png, ff, res=300, height=4, width=8, units="in")

gs_name <- 'hallmark'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_corr[[gs_name]][[1]][[4]], vasari_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), labels=c('a','b'), label_size=16, scale=.95)
dev.print(png, ff, res=300, height=6, width=10, units="in")

gs_name <- 'puchalski'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_corr[[gs_name]][[1]][[4]], vasari_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), labels=c('a','b'), label_size=16, scale=.95)
dev.print(png, ff, res=300, height=4, width=10.5, units="in")

gs_name <- 'reactome'; ff <- file.path(fdd, paste('corr', label, gs_name, 'heat.png', sep='_')); print(ff)
plot_grid(vasari_corr[[gs_name]][[1]][[4]], vasari_corr[[gs_name]][[2]][[4]], nrow=1, rel_widths=c(1, .8), labels=c('a','b'), label_size=16, scale=.95)
dev.print(png, ff, res=300, height=15, width=20, units="in")

# figure publish: vasari hallmark perts combined (gene set, single enrich, corr enrich) ----
pdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular/nn_plot_data'
s <- readRDS(file.path(pdd, 'masking_scores_vasari_hallmark.rds'))


# hallmark gene set masking, single gene masking gsea, correlatin gsea
p1 <- plot_grid(s[[1]][[4]], s[[2]][[4]], 
                NULL, NULL, NULL, NULL,
                nrow=3, rel_widths=c(.8,1), labels=c('a'), label_size=16)
p2 <- plot_grid(vasari_single[['hallmark']][[1]][[4]], vasari_single[['hallmark']][[2]][[4]], nrow=1, rel_widths=c(1, .8), labels=c('b'), label_size=16)
p3 <- plot_grid(vasari_corr[['hallmark']][[1]][[4]], vasari_corr[['hallmark']][[2]][[4]], nrow=1, rel_widths=c(1, .8), labels=c('c'), label_size=16)

ggdraw() + 
  draw_plot(p1, 0, -00, 1, 1) +
  draw_plot(p2, 0, .395, 1, .175) +
  draw_plot(p3, 0, 0, 1, .40) 
ff <- file.path(fdd, paste('combined_vasari_hallmark.png')); print(ff)
dev.print(png, ff, res=300, height=14, width=10, units="in")

# hallmark gene set masking, single gene masking gsea
p1 <- plot_grid(s[[1]][[4]], s[[2]][[4]], 
                NULL, NULL,
                nrow=2, rel_widths=c(.8,1), labels=c('a'), label_size=16)
p2 <- plot_grid(vasari_single[['hallmark']][[1]][[4]], vasari_single[['hallmark']][[2]][[4]], nrow=1, rel_widths=c(1, .8), labels=c('b'), label_size=16)

ggdraw() + 
  draw_plot(p1, 0, -.1, 1, 1.1) +
  draw_plot(p2, 0, .0, 1, .43)
ff <- file.path(fdd, paste('combined_vasari_hallmark.png')); print(ff)
dev.print(png, ff, res=300, height=10, width=12, units="in")

#
# figure 3 publish: subtype model-------
# merge single gsea and corr gsea for subtype model + coverage and single heatmap
pa <- readRDS(file.path(getwd(), 'plot_data', 'plot_subtype_single_masking.rds')) # load single gene masking heatmap
pb <- readRDS(file.path(getwd(), 'plot_data', 'plot_subtype_verhaak_gene_set_coverage.rds')) # load single gene masking coverage scatter plot

p1 <- plot_grid(s[[1]][[4]], s[[2]][[4]], nrow=1, rel_widths=c(1, .8), labels=c('a'), label_size=16)
p2 <- plot_grid(pb, nrow=1, labels=c('b'), label_size=16)

p3 <- plot_grid(subtype_single[['verhaak']][[1]][[4]], subtype_single[['verhaak']][[2]][[4]], 
                nrow=1, rel_widths=c(1, .7), scale=1, labels=c('c'), label_size=16)

p4 <- plot_grid(subtype_corr[['verhaak']][[1]][[4]], subtype_corr[['verhaak']][[2]][[4]],
                nrow=1, rel_widths=c(1, .7), scale=1, labels=c('d'), label_size=16)

ggdraw() + 
  draw_plot(p1, 0, 0, 1,1) +
  draw_plot(p2, 0.45, .7, .5, .3) +
  draw_plot(p3, 0.36, 0.35, 0.65, 0.33) +
  draw_plot(p4, 0.33, 0.0, 0.68, 0.33) 
ff <- file.path(fdd, paste('combined_subtype_single.png')); print(ff)
dev.print(png, ff, res=300, height=10.5, width=10, units="in")
