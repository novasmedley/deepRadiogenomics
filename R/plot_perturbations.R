#' After retraining, neural nets are probed with masked inputs defined by gene sets 
#' and their perturbations, or outputs are recorded in a csv file that is read in here.
#' 
#' Visualize the pertubrations
library(devtools)
library(dplyr)
library(reshape2)
library(plyr) # mapvalues
library(ggplot2)
library(gplots)
library(egg) # arrange graphs
library(pheatmap)
library(ggpubr) # title
library(cowplot) # arrange graphs
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(awtools)
library(wesanderson)
library(grid)
library(gridExtra)
library(fgsea)
library(foreach)
library(doParallel)

dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga'
fdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_figs/nn_gene_masking_r'
pdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular/nn_plot_data'
gm_dir <- file.path(dd, 'with_ae_paper_tabular','nn_gene_masking')
 
# -- functions --------------

clean_gene_ids <- function(ids) {
  ids <- gsub("-", ".", ids)
  ids <- gsub("@", ".", ids)
}

# read scores from csv files
read_gene_masking_scores <- function(gm_dir, label, type='set') {
  f <- list.files(path=gm_dir, pattern=label, full.names=F)
  f <- f[grepl('*scores',f)]
  s <- lapply(f, function(i) {
    cat('reading ', i, '\n')
    d <- read.csv(file.path(gm_dir, i), stringsAsFactors=F)
    d$name <- strsplit(i, '_')[[1]][2]
    if (d$name[1]=='canonical') {
      r <- c('KEGG','BIOCARTA','REACTOME')
      remove <- unlist(lapply(d$X, function(i) any(lapply(r, function(j) grepl(j, i))==T) ))
      d <- d[!(remove),]
    }
    d
  })
  s <- bind_rows(s)
} 
# get subset of scores for multiclass
get_values_multi <- function(scores, metric, colnames){
  a <- scores[scores$'X.1'==metric, ]
  a <- a[,c(colnames, 'X', 'name')]
  rownames(a) <- a$X
  return(a)
}
# get subset of scores for binary classes
get_values_bin <- function(scores, metric, label){
  a <- scores[, c(metric, 'X', 'name') ]
  colnames(a) <- c(label, 'X', 'name')
  rownames(a) <- a$X
  return(a)
}

# get top n genes or sets in gene masking
breaks <- seq(0, 1.0, 0.1) # legend
cols_pr <- colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(length(breaks))
cols_roc <- colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(length(breaks))

get_top_n_label <- function(n, a, a_info, type='among_all', label=NA, threshold=NA) {
  a <- as.matrix(a)
  print(n)
  print(type)
  if (type == 'among_all') {
    if (is.na(label)){
      # top n among all for all labels
      if (is.na(threshold)) thres <- sort(a, decreasing=T)[min(n, length(a))] ; cat('threshold: ', thres)
      a_thres <- a[ apply(a, 1, function(r) any(r >= thres)), , drop=F]
    } else {
      # top n for a single label
      thres <- sort(a[, label], decreasing=T)[min(n, length(a))] ; cat('threshold: ', thres)
      a_thres <- a[ a[,l] >= thres, , drop=F]
    }
  } else if (type=='among_each') {
    # among each
    # get top n gene for each label and concat
    top <- lapply(colnames(a), function(i) {
      thres <- sort(a[, i], decreasing=T)[min(n, length(a[, i]))] 
      a_thres <- a[(a[, i] >= thres), , drop=F]
      rownames(a_thres)
    })
    top <- unique(unlist(top, recursive=T))
    a_thres <- a[rownames(a) %in% top, , drop=F]
  } else {
    # among a specific label
    thres <- sort(a[, type], decreasing=T)[min(n, length(a[, type]))] 
    a_thres <- a[(a[, type] >= thres), , drop=F]
    rownames(a_thres)
  }
  
  c <- a_info[rownames(a_info) %in% rownames(a_thres), ]
  
  cat('n genes/sets', nrow(a_thres))
  return(list(mat=as.matrix(a_thres), info=c))
}

# heatmap
get_heatmaps <- function(metric, gene_set_name, pr, roc, gs_info, plot_col_names, gsname, type='among_each', top_n=20,
                         print=F, h=12, w=10, pgratio=c(1,1), fs=9, both=T, get_random=F, treeheight_row=10,
                         clus_col=T, clus_row=T, interval=100, tags=NA, show_rows=T, gs_types=c('GO'), title=NA){
  print(gene_set_name[1])
  
  if (metric=='pr_auc') {
    b <- pr
    b2 <- roc
    cols1 <- cols_pr
    cols2 <- cols_roc
    met_name1 <- 'average precision'
    met_name2 <- 'AUC'
    ctitle1 <- 'average precision'
    ctitle2 <- 'AUC'
  } else {
    b <- roc
    b2 <- pr
    cols1 <- cols_roc
    cols2 <- cols_pr
    met_name1 <- 'AUC'
    met_name2 <- 'average precision'
    ctitle1 <- ''
    ctitle2 <- 'average precision'
  }
  met_name <- 'score\nrange'
  
  colornames <- as.character(round(seq(0,1,.1), digits=1))
  names(cols1) <- colornames
  names(cols2) <- colornames
  
  # keep certain gene sets
  if (length(gene_set_name) > 1) {
    b <- b[b$name %in% gs_types, ]
    b2 <- b2[b2$name %in% gs_types, ]
  }
  
  # get row info
  info <- b[, c('X','name'), drop=F]
  info[,'gene set'] <- info$name
  
  info <- left_join(info, verhaak_labels[, c('X','subtype')]) # attach gene set data
  info$subtype <- 'na'
  info$subtype[info$X %in% verhaak_labels$X[verhaak_labels$subtype=='PN']] <- 'PN'
  info$subtype[info$X %in% verhaak_labels$X[verhaak_labels$subtype=='CL']] <- 'CL'
  info$subtype[info$X %in% verhaak_labels$X[verhaak_labels$subtype=='MES']] <- 'MES'
  info$subtype[info$X %in% verhaak_labels$X[verhaak_labels$subtype=='NL']] <- 'NL'
  info$subtype[info$X %in% verhaak_labels$X[verhaak_labels$subtype=='']] <- 'unnamed'
  
  info <- left_join(info, gs_info)
  rownames(info) <- rownames(b) # attach back rownames, must have left_joins with info on left always
  
  b <- b[, plot_col_names, drop=F]
  b2 <- b2[, plot_col_names, drop=F]
  
  if (length(gene_set_name) > 1) {
    a <- b[info$X %in% gene_set_name, , drop=F] # for searched gs, for a specific gene set
    cat('gs left: ', nrow(a))
  } else {
    # for a gs collection
    gsname <- gene_set_name
    if (gene_set_name=='all') {
      a <- b # check top gene sets for all types of gene sets
    } else {
      if (gene_set_name=='verhaak' & get_random==T) {
        a <- b[info$name %in% c('verhaak','random'), ] # get random too
      } else {
        a <- b[info$name==gene_set_name, ] # for a specific gene set type
      }
    }
  }
  

  pdat <- get_top_n_label(n=top_n, a=a, a_info=info, type=type)
  if (!is.na(tags)){
    rn <- rownames(pdat$info)
    pdat$info <- left_join(pdat$info, tags, copy=T)
    rownames(pdat$info) <- rn # attach back names
  }
  pdat2 <- pdat
  pdat2$mat <- as.matrix(b2[rownames(pdat$mat), colnames(pdat$mat), drop=F]) # get corresponding data in other metric 
  
  
  if (nrow(pdat$mat)==1) clus_row <- F
  if (ncol(pdat$mat)==1) clus_col <- F
  
  if (gene_set_name[1]=='verhaak') {interval <- 50}
  if (gene_set_name[1]=='puchalski') {
    if ('Neural' %in% plot_col_names) {
      # match naming in paper
      f <- pdat$info
      m <- pdat$mat
      f$X <- gsub('^Barres', 'Zhang', f$X)
      f$X <- gsub('^DarmanisBarres', 'Darmanis', f$X)
      rownames(f) <- f$X
      rownames(m) <- f$X
      # match order of rows and cols as in their paper
      paper_rows <- c('Zhang_matureastrocytes', 'Zhang_neuron', 'Zhang_oligodendrocytes', 
                      'Darmanis_astrocytes', 'Darmanis_mixOPCOligNeurons','Zhang_endothelial',
                      'Patel_immune','Zhang_microgliamacrophages','Patel_anticellcycle','Patel_hypoxia')
      paper_cols <- c('Neural', 'Proneural','Classical', 'Mesenchymal')
      
      m <- rbind(m[paper_rows, paper_cols], m[!(rownames(m) %in% paper_rows), paper_cols])
      f <- f[rownames(m),] # match in info too
      pdat$info <- f
      pdat$mat <- m
      
      m2 <- pdat2$mat
      rownames(m2) <- f$X
      m2 <- rbind(m2[paper_rows , paper_cols], m2[!(rownames(m) %in% paper_rows), paper_cols])
      pdat2$mat <- m2
      pdat2$info <- f
      
      clus_col <- F; clus_row <- F
    } else {
      # match naming in paper
      f <- pdat$info
      m <- pdat$mat
      f$X <- gsub('^Barres', 'Zhang', f$X)
      f$X <- gsub('^DarmanisBarres', 'Darmanis', f$X)
      rownames(f) <- f$X
      rownames(m) <- f$X
      pdat$info <- f
      pdat$mat <- m
      pdat2$info <- f
      pdat2$mat <- m
    }
  }

  
  if (gene_set_name[1]=='single') {
    ann_row <- pdat$info[, 'subtype', drop=F]
    ann_colors <- list('subtype'=subtype_color)
    pt <- 'genes'
  } else {
    # info$cover_new <- as.character(round_any(info$cover, 10, f=floor)) # round to tens
    # info$size_new <- as.character(round_any(info$size, 100, f=floor)) # round to hundreds
    pdat$info$cover_new <- cut(pdat$info$cover,
                               breaks=seq(0,100,10)
                               # labels=paste('<', seq(10,100,10))
    )
    
    maxsize <- max(na.omit(pdat$info$size))
    if (maxsize <= 100) interval <- 10
    pdat$info$size_new <- cut(pdat$info$size, 
                              breaks=seq(0,maxsize+interval, interval)
                              # labels=paste('<', seq(100,maxsize,interval))
    )
    pdat2$info$cover_new <- pdat$info$cover_new
    pdat2$info$size_new <- pdat$info$size_new
    
    nc <- levels(droplevels(pdat$info$cover_new)) 
    gene_cover_color <- rev(colorRampPalette(rev(brewer.pal(n=6, name='Greens')))(length(nc)))
    names(gene_cover_color) <- nc
    
    ns <- levels(droplevels(pdat$info$size_new)) 
    gene_size_color <- colorRampPalette(brewer.pal(n=6, name="Blues" ))(length(ns))
    names(gene_size_color) <- ns
    
    pdat$info$'% coverage' <- as.character(pdat$info$cover_new)
    pdat$info$'num. genes' <- as.character(pdat$info$size_new)
    
    if (!(is.na(tags))) {
      gt <- unique(pdat$info$query)
      ann_colnames <- c('gene set','num. genes','% coverage', 'query')
      ann_row <- pdat$info[, rev(ann_colnames)]
      ann_colors <- list('gene set'=gene_set_color[names(gene_set_color) %in% pdat$info$'gene set'],
                         '% coverage'=gene_cover_color,
                         'num. genes'=gene_size_color,
                         'query'=gene_tag_color[names(gene_tag_color) %in% pdat$info$'query'])
      
    } else {
      ann_colnames <- c('gene set', 'num. genes','% coverage')
      ann_row <- pdat$info[, rev(ann_colnames)]
      ann_colors <- list('gene set'=gene_set_color[names(gene_set_color) %in% pdat$info$'gene set'],
                         '% coverage'=gene_cover_color,
                         'num. genes'=gene_size_color)
    }

    

    pt <- 'gene sets'
  }
  
  # rename for plotting
  rn <- gene_set_name[1] %in% c('verhaak', 'single', 'puchalski') # rename for plotting
  lr <- rownames(pdat$info)
  lr[!(rn)] <- tolower(lr[!(rn)]) 
  lr <- gsub('hallmark_','', lr)
  lr <- gsub('reactome_','', lr)
  lr <- gsub('go_','', lr)
  lr <- gsub('biocarta_','', lr)
  lr <- gsub('kegg_','', lr)
  lr <- gsub('_',' ', lr)
  
  rn <- pdat$info$'name' %in% c('immuno') # remove first word
  if (any(rn==T)) {
    lr[rn] <- unlist(lapply(lr[rn], function(i) {
      k <- strsplit(i, ' ')[[1]]
      k <- paste(k[2:length(k)], collapse=' ')
    }))
  }
  
  lr <- lapply(lr, function(i) { # shorten long names
    k <- strsplit(i, ' ')[[1]]
    # if (nchar(i) > 40)  {
    #   k <- paste(k[1:5], collapse=' ') # truncate to first 5 words
    #   k <- paste0(k, '...')
    if (nchar(i) > 40 & nchar(i) < 80)  {
    p <- floor(length(k)/2)
      k <- paste0(paste(k[1:p], collapse=' '), '\n', paste(k[(p+1):length(k)], collapse=' '))
    } else if (nchar(i) >= 80) {
      p <- floor(length(k)/3)
      k <-  k <- paste0(paste(k[1:p], collapse=' '), '\n', paste(k[(p+1):(p+p)], collapse=' ')) # truncate
    } else{
      k <- i
    }
    k
  })
  lr <- unlist(lr)
  
  hp1 <- pheatmap(mat=pdat$mat,
                  main=met_name1,
                  annotation_row=ann_row, annotation_colors=ann_colors, annotation_legend=F,
                  cluster_cols=clus_col, cluster_rows=clus_row, fontsize=fs, fontsize_row=fs+2, fontsize_col=fs+2,
                  legend=F, lineheight=.7,
                  color=cols1, breaks=breaks, labels_row=lr, 
                  show_colnames=T, show_rownames=show_rows, border_color = "grey60",
                  treeheight_col=0, treeheight_row=treeheight_row)
  if ((gene_set_name[1]=='puchalski' & 'Neural' %in% plot_col_names)) {
    mat2 <- pdat2$mat
    lr2 <- lr
  } else {
    lr2 <- lr[hp1$tree_row$order] # rename according to reorder
    if (clus_col==F) {
      mat2 <- pdat2$mat[hp1$tree_row$order, , drop=F] # reorder values based on clustering in hp1
    } else {
      mat2 <- pdat2$mat[hp1$tree_row$order, hp1$tree_col$order] # reorder values based on clustering in hp1
    }
  }
  hp2 <- pheatmap(mat=mat2, 
                  main=met_name2,
                  annotation_row=ann_row, annotation_colors=ann_colors, annotation_legend=T,
                  cluster_cols=F, cluster_rows=F,  fontsize=fs, fontsize_row=fs+2, fontsize_col=fs+2,
                  legend_labels=breaks, legend=T, legend_breaks=breaks, drop_levels=T, lineheight=.7,
                  color=cols2, breaks=breaks, labels_row=lr2, 
                  show_colnames=T, show_rownames=show_rows, border_color = "grey60")
  
  if (both) {
    fname <- paste(lt, gsname ,type,'top', as.character(top_n), metric, as.character(nrow(pdat$mat)), sep='_'); print(fname)
    
    dev.new();p <- plot_grid(hp1[[4]], hp2[[4]], rel_widths=pgratio); print(p)
    if (print) dev.print(png, file.path(fdd, paste0(fname,'_both','.png')), res=300, height=h, width=w, units="in")
    
    l <- list(hp1, hp2)
    names(l) <- c(met_name1, met_name2)
  } else {
    fname <- paste(lt, gsname ,type,'top', as.character(top_n), metric, as.character(nrow(pdat$mat)), sep='_'); print(fname)
    
    p <- pheatmap(mat=pdat$mat,
                  main=title,
                  annotation_row=ann_row, annotation_colors=ann_colors, annotation_legend=T,
                  cluster_cols=clus_col, cluster_rows=clus_row, fontsize=fs, fontsize_row=fs+2, fontsize_col=fs+2,
                  legend=T, lineheight=.7,
                  color=cols1, breaks=breaks, labels_row=lr, 
                  show_colnames=T, show_rownames=show_rows, border_color = "grey60",
                  treeheight_col=0, treeheight_row=treeheight_row); dev.new();print(p)
    if (print) dev.print(png, file.path(fdd, paste0(fname,'.png')), res=300, height=h, width=w/1.5, units="in")
    l <- list(p)
    names(l) <- c(met_name1)
  }
  
  l
}
#
# -- patient info  --------------

cdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'
clin <- 'nationwidechildrens.org_clinical_patient_gbm.txt'
clin <- read.table(file.path(cdd, clin), sep='\t', skip=1, header=T, stringsAsFactors=F) # get clinical
clin <- clin[-1,]
clin$X <- clin$bcr_patient_barcode
clin$age_cat <- factor(unlist(lapply(clin$age_at_initial_pathologic_diagnosis, function(i) ifelse(i < 50, 1, 0))))
clin$age_at_initial_pathologic_diagnosis <- as.numeric(clin$age_at_initial_pathologic_diagnosis)
clin$karnofsky_performance_score <- as.numeric(clin$karnofsky_performance_score)
clin$eastern_cancer_oncology_group <- as.numeric(clin$eastern_cancer_oncology_group)

# load patient label info
vasari <- 'vasari_annotations_191.csv'; vasari <- read.csv(file.path(cdd, vasari), row.names=1)
vasari$X <- rownames(vasari)


# -- gene set info  --------------

gs_dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'

verhaak <- read.table(file.path(gs_dd, 'TCGA_unified_CORE_ClaNC840.txt'), sep='\t', header=T, stringsAsFactors=F)
verhaak <- verhaak[-1,1:2]
colnames(verhaak) <- c('go_genes', 'subtype')
verhaak_labels <- verhaak
colnames(verhaak_labels) <- c('X', 'subtype')

verhaak$X <- verhaak$go_genes
t <- lapply(c('PN','NL','CL','MES'), function(i) verhaak$go_genes[verhaak$subtype==i])
t[['all']] <- verhaak$go_genes
verhaak <- t
names(verhaak) <- c('PN','NL','CL','MES', 'all')

puchalski <- read.table(file.path(gs_dd, 'gene_sets_Puchalski', 'aaf2666_Table-S15.csv'), sep=',', header=T, stringsAsFactors=F)
puchalski <- puchalski[-1,]
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
gids <- clean_gene_ids(read.csv(file.path(dd,'with_ae_paper_tabular','nn_gene_masking', 'subtype_single_roc_scores.csv'), stringsAsFactors=F)$X)

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

rand_info <- data.frame(cover=c(100, 100), size=c(100,200), X=c('random_100','random_200'))
rand_info$'gene set' <- c('random','random')
gs_info <- bind_rows(gs_info, rand_info)
rm(rand_info)

# # remove duplicates
# gs_info <- split(gs_info, gs_info$`gene set`)
# gs_info$canonical <- gs_info$canonical[ !(gs_info$canonical$X %in% c(gs_info$biocarta$X,gs_info$kegg$X, gs_info$reactome$X)), ]
# gs_info <- do.call('rbind', gs_info)
# rownames(gs_info) <- gs_info$X

# -- heatmap colors, global -----
subtype_color <- c('white','grey','red','yellow','blue','green')
names(subtype_color) <- c('na','unnamed','NL','PN','MES','CL')

gene_set_color <- c(RColorBrewer::brewer.pal(name='Paired',n=12), RColorBrewer::brewer.pal(name='Set3',n=4))
names(gene_set_color) <- c('hallmark','chromosome','motif','computational','GO',
                           'onco','immuno', 'verhaak', 'puchalski', 'single',
                           'reactome', 'random',
                           'kegg','biocarta','canonical','chem')
unique(gs_info$`gene set`) %in% names(gene_set_color)



#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# Verhaak model + vasari gene set pertubrations  ----------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------

# veer <- read.csv(file.path(dd,'experiments/gbm_veerhak_retrain_18_12_05/veerhak_perturbations.csv'))
vdd <- file.path(dd,'with_ae_paper_tabular/nn_gene_masking')
vpreds <- c('subtype_random_verh_random_200_preds.csv', 'subtype_random_verh_random_100_preds.csv',
            'subtype_verhaak_all_preds.csv', 
            'subtype_verhaak_PN_preds.csv','subtype_verhaak_NL_preds.csv',
            'subtype_verhaak_MES_preds.csv', 'subtype_verhaak_CL_preds.csv')
names(vpreds) <- c('random_200','random_100','all','PN','NL','MES','CL')
veer <- lapply(names(vpreds), function(i) {
  x <- read.csv(file.path(vdd, vpreds[[i]])) 
  x$gene_set <- i
  x
})
veer <- bind_rows(veer)
colnames(veer)
veer <- veer[, !(colnames(veer) %in% c('y_true'))]
veer <- melt(veer)
colnames(veer) <- c('id','subtype','gene_set','class_prob','prediction')

# get sorted heatmap
order_by = 'MES.probability'
veer <- split(veer, f=veer$gene_set)

cat(names(veer))  # order_by should match names below
order_by <- c('MES.probability', 'CL.probability', 'MES.probability', 'NL.probability', 'PN.probability', 'MES.probability')

i <- 0
veer <- lapply(veer, function(gs){
  df <- split(gs, f=gs$subtype)
  i <<- i + 1
  this <- order_by[i]
  cat(this)
  df <- lapply(df, function(s){  # get id order
    s <- s[s$class_prob==this,]
    s <- s[order(s$prediction),]
    s$order <- rownames(s)
    s$id
  })
  
  gs$id <- factor(gs$id, levels=unlist(df)) # set id order as level order
  gs$class_prob <- gsub('.probability', '', gs$class_prob)
  gs
})

# must plot each df in veer separately due to diff id level ordering
dev.new()

#' plot contants
fs <- 9
pTheme <- theme(legend.position='none',
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.title.y=element_text(margin=margin(0,-0.2,0,0), angle=0, size=fs, color="#666666"),
                axis.title.x=element_blank(),
                axis.ticks=element_blank(),
                axis.line=element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.margin=unit(c(0,0,0,0), "mm") )
x <- scale_x_discrete(expand=c(0,0))
y <- scale_y_discrete(expand=c(0,0))
fill_s <- scale_fill_viridis(limits=c(0,1),breaks=seq(0, 1, 0.2),begin=0, end=1)

plot_veerhak_perturbs <- function(gs, yname=names(gs)) {
  gs <- split(gs, gs$subtype)
  cat(names(gs))
  names(gs) <- c('CL','MES','NL','PN')
  cat(names(gs))
  
  p4 <- ggplot(gs[[4]], aes(x=class_prob, y=id, fill=prediction)) + geom_tile() + pTheme + x + y + fill_s +
    ylab(yname[4]) +
    theme(axis.text.x=element_text(size=fs, color="#666666"))
  
  p3 <- ggplot(gs[[3]], aes(x=class_prob, y=id, fill=prediction)) + geom_tile() + pTheme + x + y + fill_s+
    ylab(yname[3])
  
  p2 <- ggplot(gs[[2]], aes(x=class_prob, y=id, fill=prediction)) + geom_tile() + pTheme + x + y + fill_s+
    ylab(yname[2])
  
  p1 <- ggplot(gs[[1]], aes(x=class_prob, y=id, fill=prediction)) + geom_tile() + pTheme + x + y + fill_s +
    ylab(yname[1])

  a <- egg::ggarrange(p1, p2, p3, p4, nrow=4)
}




# GET legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  } 
this <- 2
names(veer)[this]
gs <- veer[[this]]
gs <- split(gs, gs$subtype)
p3 <- ggplot(gs[[2]], aes(x=class_prob, y=id, fill=prediction)) + geom_tile() + pTheme + x + y + fill_s+
  ylab(names(gs)[2]) + 
  theme(legend.position='right', 
        legend.title=element_text(size=fs),
        legend.text=element_text(size=fs, color="grey30"),
        legend.margin=margin(l = .15, unit='cm'),
        legend.key.size = unit(.7, "cm")) + 
  labs(fill="probability")
veer_legend <- g_legend(p3)




# plot heatmaps
veer_plots <- lapply(veer, function(i) plot_veerhak_perturbs(i))
veer_plots <- lapply(veer_plots, function(i) annotate_figure(i, top = text_grob(" ", size = fs+2), left = text_grob("truth", size = fs+2)))

veer_plots2 <- lapply(veer[1:5], function(i) plot_veerhak_perturbs(i, yname=rep('',4) ))
veer_plots2 <- lapply(veer_plots2[1:5], function(i) annotate_figure(i, top = text_grob(" ", size = fs+2)))


names(veer_plots) # order of naming
ps <- plot_grid(veer_plots[[6]], veer_plots2[[1]], veer_plots2[[2]], veer_plots2[[3]], veer_plots2[[4]], veer_plots2[[5]], veer_legend, 
                labels=c("\t\t   random 200 genes", " all subtype genes", " CL genes", " MES genes", " NL genes", " PN genes"),
                # labels=c("a", "b", "c", "d", "e", "f"),
                nrow=1, rel_widths=c(1.6,1,1,1,1,1,.5),
                hjust=0, vjust=1.5,
                label_size=fs)
ps <- annotate_figure(ps, bottom=text_grob(' model prediction',size=fs+2))
print(ps)
dev.print(png, file.path(fdd, paste0('subtype_self_set_masking','.png')), res=300, height=6, width=13, units="in")

ps_few <- plot_grid(veer_plots[[6]], veer_plots2[[1]], veer_plots2[[3]], veer_legend, 
                labels=c("\t\t       random 200 genes", " all subtype genes", " MES genes"),
                # labels=c("a", "b", "c", "d", "e", "f"),
                nrow=1, rel_widths=c(1.5,1,1,.5),
                hjust=0, vjust=1.5,
                label_size=fs)
ps_few <- annotate_figure(ps_few, bottom=text_grob(' model prediction',size=fs+2))
print(ps_few)
dev.print(png, file.path(fdd, paste0('subtype_self_set_masking_fewer','.png')), res=300, height=6, width=9, units="in")

#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------


# Verhaak subtype single gene set masking plot coverage of top n genes in Verhaak  -------
sub_nickname <- c('CL','MES','NL','PN')
sub_names <- c('Classical','Mesenchymal','Neural','Proneural')
lt <- 'subtype'

metric <- 'pr_auc'
mname <- 'AP'
subtype <- get_values_multi(scores=read_gene_masking_scores(gm_dir, label='subtype'), metric=metric, colnames=sub_names)

a <- as.matrix(subtype[subtype$name=='single', sub_names])
length(rownames(a)[rownames(a) %in% verhaak$all])/length(verhaak$all) # ratio of verhaak genes in microarray



parallel::mcaffinity(1:7)
cl <- makeCluster(7)
registerDoParallel(cl)
n <- 12042
n_in_verhaak <- foreach (i=1:n) %dopar% {
  thres <- sort(a, decreasing=T)[i] 
  a2 <- a[apply(a, 1, function(r) any(r >= thres)), ]
  
  topg <- rownames(a2) # check if top genes are in verhaak gene set
  length(topg[topg %in% verhaak$all]) # number of top n genes in verhaak genes
  length(verhaak$all)
  x1 <- 100*length(topg[topg %in% verhaak$all])/length(verhaak$all)
  x2 <- 100*length(topg[topg %in% verhaak$all])/i
  data.frame(in_verh=x1, in_top=x2)
}
stopCluster(cl)

t <- bind_rows(n_in_verhaak)
t <- melt(t)
t$x <- seq(1,n,1)/1000
saveRDS(t, file.path(getwd(), 'plot_data', 'subtype_verhaak_gene_set_coverage.rds'))

t <- readRDS(file.path(getwd(), 'plot_data', 'subtype_verhaak_gene_set_coverage.rds'))
topn_sub <- ggplot(t[t$variable %in% c('in_verh','in_top') ,], aes(x=x, y=value, color=variable)) + 
  geom_point(size=0.1) +
  scale_color_manual(labels = c("% subtype genes covered", "% top N that are subtype genes"), values = c("blue", "red"))  +
  scale_y_continuous(breaks=seq(0,100,20)) +
  scale_x_continuous(breaks=seq(0,13, 1), minor_breaks=seq(0.5,12.5,.5)) +
  labs(x=bquote('top N genes ranked by average precision '~(10^3)), y='%') +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme_bw() +  theme(legend.justification = c(1, 1), 
                      legend.position = c(0.95, 0.75),
                      legend.title=element_blank(),
                      legend.text=element_text(size=fs),
                      axis.text=element_text(size=fs),
                      axis.title=element_text(size=fs, angle=0)); print(topn_sub)
saveRDS(topn_sub, file.path(getwd(), 'plot_data', 'plot_subtype_verhaak_gene_set_coverage.rds'))
dev.print(png, file.path(fdd, paste0('subtype_',gs_type, '_masking_', metric_name,'_curves','.png')), res=300, height=2, width=8, units="in")






# Verhaak subtype gene set masking scores ---------------------------------------------------------------------------

# pheatmap

sub_nickname <- c('CL','MES','NL','PN')
sub_names <- c('Classical','Mesenchymal','Neural','Proneural')
lt <- 'subtype'

subtype_pr <- get_values_multi(scores=read_gene_masking_scores(gm_dir, label='subtype'), metric='pr_auc', colnames=sub_names)
subtype_roc <- get_values_multi(scores=read_gene_masking_scores(gm_dir, label='subtype'), metric='roc_auc', colnames=sub_names)

unique(subtype_pr$name)

# invididual gs plots

lt <- 'subtype'
b <- T
s_verh <- get_heatmaps(metric='pr_auc', gene_set_name='verhaak', pr=subtype_pr, roc=subtype_roc, gs_info=gs_info, plot_col_names=sub_names, type='among_each', top_n=20,
                  print=T, h=2.5, w=6, pgratio=c(.7,1), fs=9, both=F, get_random=T)
s <- get_heatmaps(metric='pr_auc', gene_set_name='puchalski', pr=subtype_pr, roc=subtype_roc, gs_info=gs_info, plot_col_names=sub_names, type='among_each',
                  print=T, h=5, w=10, pgratio=c(.8,1), fs=8, both=b, top_n=50) # take all
s <- get_heatmaps(metric='pr_auc', gene_set_name='hallmark', pr=subtype_pr, roc=subtype_roc, gs_info=gs_info, plot_col_names=sub_names, type='among_each',
                  print=T, h=6, w=10, pgratio=c(.8,1), fs=8, both=b, top_n=20) 
s_single <- get_heatmaps(metric='pr_auc', gene_set_name='single', pr=subtype_pr, roc=subtype_roc, gs_info=gs_info, plot_col_names=sub_names, type='among_each', top_n=20,
                  print=T, h=12, w=8, pgratio=c(.75,1), fs=8, both=F)
s <- get_heatmaps(metric='pr_auc', gene_set_name='chromosome', pr=subtype_pr, roc=subtype_roc, gs_info=gs_info, plot_col_names=sub_names, type='among_each', top_n=20,
                  print=T, h=12, w=8, pgratio=c(.75,1), fs=11, both=b) 
s <- get_heatmaps(metric='pr_auc', gene_set_name='onco', pr=subtype_pr, roc=subtype_roc, gs_info=gs_info, plot_col_names=sub_names, type='among_each', top_n=20,
                  print=T, h=8, w=10, pgratio=c(.8,1), fs=11, both=b) 

# save single for combinging figs in plot_gsea.R
saveRDS(s_single, file.path(getwd(), 'plot_data', 'plot_subtype_single_masking.rds'))

# combined gs plot
s <- get_heatmaps(metric='pr_auc', gene_set_name='all', pr=subtype_pr, roc=subtype_roc, gs_info=gs_info, plot_col_names=sub_names,
                  type='among_all', top_n=5,
                  print=T, # not helpful
                  h=30, w=15, pgratio=c(.4,1), fs=11, both=b) 


# combine subtype plots, Fig 2 subtype masking  ---------------------------------------------------------------------------
dev.new();
p1 <- plot_grid(ps_few, NULL, nrow=2, rel_heights=c(1,.6),labels=c('a'), label_size=16)
p2 <- plot_grid(NULL,s_verh[[1]][[4]], nrow=1, rel_widths=c(0.1,1),labels=c('b'), label_size=16)
ggdraw() + 
  draw_plot(p1, 0,0, 1,1) +
  draw_plot(p2, 0, 0, .93,.37)
dev.print(png, file.path(fdd, paste0('subtype_masking_combined_v','.png')), res=300, height=8, width=6, units="in")

## TODO get the GO annotations Verhaak did

# get top n genes for each label and combine into one  ----
b <- subtype
lt <- 'subtype'

info <- b[, c('X','name'), drop=F]
info[,'gene set'] <- info$name

info <- left_join(info, verhaak[, c('X','subtype')])

singlerows <- info$x[info$name=='single']
names(subtype_color)
info$subtype <- 'na'
info$subtype[info$X %in% verhaak$go_genes[verhaak$subtype=='PN']] <- 'PN'
info$subtype[info$X %in% verhaak$go_genes[verhaak$subtype=='CL']] <- 'CL'
info$subtype[info$X %in% verhaak$go_genes[verhaak$subtype=='MES']] <- 'MES'
info$subtype[info$X %in% verhaak$go_genes[verhaak$subtype=='NL']] <- 'NL'
info$subtype[info$X %in% verhaak$go_genes[verhaak$subtype=='']] <- 'unnamed'

b <- b[, sub_names]

rownames(b) <- tolower(rownames(b)) # rename for plotting
rownames(b) <- gsub('hallmark_','',rownames(b))
rownames(info) <- rownames(b)
n <- 50; fname <- paste(lt, 'among_each_top', as.character(n), metric, sep='_'); fname;
pdat <- get_top_n_label(n=n, a=b[!(info$name=='single'), ], a_info=info, type='among_each')
pheatmap(mat=pdat$mat, main=paste(metric, as.character(nrow(pdat$mat)), 'rows'),
         cluster_cols=T, color=cols, breaks=breaks, show_colnames=T, show_rownames=T, legend=T, border_color = "grey60",
         annotation_row=pdat$info[, 'gene set', drop=F], annotation_colors=list('gene set'=gene_set_color), treeheight_col=20, treeheight_row=20)
dev.print(png, file.path(fdd, paste0(fname,'.png')), res=300, height=12, width=10, units="in")

# get top n genes among all labels for each gene set



#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
# Vasari models + all gene set masking scores  ---------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------
#' -----------------------------------------------------------------------------------------------------------------------------


# get data -------
labels <- c('f5','f6','f7','f9','f10','f14')
lt <- 'vasari'

get_vasari_gs_scores <- function(metric) {
  vasari_set <- c(lapply(labels, function(i) get_values_bin(scores=read_gene_masking_scores(gm_dir, label=i), 
                                                              metric=metric,
                                                              label=i)))
  vasari_set <- Reduce(function(x,y) left_join(x,y, by=c('X','name')), vasari_set)
  
  
  # remove large gene sets
  vasari_set <- split(vasari_set, vasari_set$name)
  vasari_set <- lapply(names(vasari_set), function(i){
    b <- vasari_set[[i]]
    if (i != 'single'){
      keep <- gs_info$X[gs_info$size < 500 & gs_info$`gene set`==i]
      b <- b[b$X %in% keep, ]
    }
    b
  })
  vasari_set <- bind_rows(vasari_set)
  rownames(vasari_set) <- vasari_set$X
  vasari_set
}
vasari_pr <- get_vasari_gs_scores(metric='pr_auc')
vasari_roc <- get_vasari_gs_scores(metric='roc_auc')

# rename for plotting
old <- c('f5', 'f6', 'f7', 'f9', 'f10', 'f14')
new <- c('enhancing', 'nCET','necrosis','focal', 'infiltrative','edema')
colnames(vasari_pr) <- mapvalues(colnames(vasari_pr),
                       from=c(old,'X'),
                       to =c(new, 'X'))
colnames(vasari_roc) <- mapvalues(colnames(vasari_roc),
                        from=c(old, 'X'),
                       to =c(new, 'X'))
labels <- c('nCET', 'enhancing','edema','focal', 'necrosis','infiltrative') # ordering if not clustering cols in heatmaps
unique(vasari_pr$name) # available gs data

# pheatmap of gene set scores  -------

# invididual gs plots
# test
metric='pr_auc'
gene_set_name='single'
pr=vasari_pr
roc=vasari_roc
gs_info=gs_info
plot_col_names=labels
type='among_each'
top_n=10
print=F
h=12
w=10
pgratio=c(1,1)
fs=9
both=T
get_random=F
treeheight_row=10
clus_col=T
clus_row=T
interval=100
tags=NA
show_rows=T
gs_types=c('GO')
# test

lt <- 'vasari'
b <- T
s <- get_heatmaps(metric='pr_auc', gene_set_name='verhaak', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=2.5, w=6, pgratio=c(.75,1), fs=9, both=b)
s <- get_heatmaps(metric='pr_auc', gene_set_name='puchalski', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=4, w=9, pgratio=c(.85,1), fs=7, both=b, top_n=20) # get all
s <- get_heatmaps(metric='pr_auc', gene_set_name='hallmark', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=5, w=8, pgratio=c(.8,1), fs=9, both=b, top_n=5)
saveRDS(s, file.path(pdd, 'masking_scores_vasari_hallmark.rds'))
s <- get_heatmaps(metric='pr_auc', gene_set_name='single', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=8.5, w=8, pgratio=c(.8,1), fs=9, both=b, top_n=10)
s <- get_heatmaps(metric='pr_auc', gene_set_name='canonical', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=3.5, w=9, pgratio=c(.8,1), fs=6, both=b, top_n=5)

s <- get_heatmaps(metric='pr_auc', gene_set_name='chromosome', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=6, w=7, pgratio=c(.75,1), fs=8, both=b, top_n=5)
s <- get_heatmaps(metric='pr_auc', gene_set_name='onco', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=6, w=10, pgratio=c(.8,1), fs=, both=b, top_n=5)
s <- get_heatmaps(metric='pr_auc', gene_set_name='GO', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=8, w=10, pgratio=c(.8,1), fs=9, both=b, top_n=5)

s <- get_heatmaps(metric='pr_auc', gene_set_name='reactome', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=8, w=10, pgratio=c(.8,1), fs=9, both=b, top_n=5)
s <- get_heatmaps(metric='pr_auc', gene_set_name='motif', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  print=T, h=6, w=10, pgratio=c(.8,1), fs=8, both=b, top_n=5)

# did not get data
# s <- get_heatmaps(metric='pr_auc', gene_set_name='motif', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
#                   print=T, h=8, w=10, pgratio=c(.8,1), fs=11, both=b, top_n=5)
# s <- get_heatmaps(metric='pr_auc', gene_set_name='immuno', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
#                   print=T, h=8, w=10, pgratio=c(.4,1), fs=11, both=b, top_n=5)


# combined gs plot
s <- get_heatmaps(metric='pr_auc', gene_set_name='all', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  type='among_each', top_n=5,
                  print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b)
s <- get_heatmaps(metric='pr_auc', gene_set_name='all', pr=vasari_pr, roc=vasari_roc, gs_info=gs_info, plot_col_names=labels,
                  type='among_all', top_n=50,
                  print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b) 

s <- get_heatmaps(metric='pr_auc', gene_set_name='all', 
                  pr=vasari_pr[!(vasari_pr$name %in% c('chem', 'onco','motif')),], roc=vasari_roc[!(vasari_roc$name %in% c('chem', 'onco','motif')),], 
                  gs_info=gs_info, plot_col_names=labels,print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b,
                  type='enhancing', top_n=15)
s <- get_heatmaps(metric='pr_auc', gene_set_name='all', 
                  pr=vasari_pr[!(vasari_pr$name %in% c('chem', 'onco','motif')),], roc=vasari_roc[!(vasari_roc$name %in% c('chem', 'onco','motif')),], 
                  gs_info=gs_info, plot_col_names=labels,print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b,
                  type='edema', top_n=15)

s <- get_heatmaps(metric='pr_auc', gene_set_name='all', 
                  pr=vasari_pr[!(vasari_pr$name %in% c('chem', 'onco','motif')),], roc=vasari_roc[!(vasari_roc$name %in% c('chem', 'onco','motif')),], 
                  gs_info=gs_info, plot_col_names=labels,print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b,
                  type='focal', top_n=20)
s <- get_heatmaps(metric='pr_auc', gene_set_name='all', 
                  pr=vasari_pr[!(vasari_pr$name %in% c('chem', 'onco','motif')),], roc=vasari_roc[!(vasari_roc$name %in% c('chem', 'onco','motif')),], 
                  gs_info=gs_info, plot_col_names=labels,print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b,
                  type='nCET', top_n=20)

# none above 0.7
s <- get_heatmaps(metric='pr_auc', gene_set_name='all', 
                  pr=vasari_pr[!(vasari_pr$name %in% c('chem', 'onco','motif')),], roc=vasari_roc[!(vasari_roc$name %in% c('chem', 'onco','motif')),], 
                  gs_info=gs_info, plot_col_names=labels,print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b,
                  type='necrosis', top_n=20)
s <- get_heatmaps(metric='pr_auc', gene_set_name='all', 
                  pr=vasari_pr[!(vasari_pr$name %in% c('chem', 'onco','motif')),], roc=vasari_roc[!(vasari_roc$name %in% c('chem', 'onco','motif')),], 
                  gs_info=gs_info, plot_col_names=labels,print=T, h=6, w=15, pgratio=c(.8,1), fs=9, both=b,
                  type='infiltrative', top_n=20)

# get literature related findings -------
## radiogenomics

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
parsons <- c('TP53', 'PI3K','RB1','EGFR','PTEN','CDKN2A', 'NF1','CDK4','IDH1','PIK3CA', 'PIK3R1') 
mclendon <- c('EGFR', 'CDK4', 'PDGFRA','MDM2','MET','CDK6','MYCN','CCND2','PIK3CA','AKT3', 'CDKN2A',
              'CDKN2B','CDKN2C','RB1','PARK2','NF1','ERBB2') 
gen <- unique(c(parsons,mclendon))

## function to read saved query results
get_search <- function(fname, tag) {
  print(fname)
  gs <- gmtPathways(gmt.file=file.path(dd,'model_ready','msigdb_v6.2_GMTs','search_results', paste0(fname,'.gmt')))
  cat('num gs: ', length(gs))
  gs <- data.frame(X=na.omit(names(gs)), query=fname)
}

labels <- c('enhancing', 'nCET','necrosis','focal', 'infiltrative','edema') # don't plot useless labels

## get lit heatmaps  -----
# enhan <- c(diehn_enhan, jam_enhan_genes, gut_enhan)
# enhan <- lapply(enhan, get_search); enhan <- bind_rows(enhan)
# enhan <- rbind(enhan, data.frame(X=jam_enhan_paths, query='named'))
# enhan <- enhan[!(duplicated(enhan$X)),] # if using as tags, need to remove duplicates
# cat('total gs: ', nrow(enhan))
# enhan_p <- get_heatmaps(metric='pr_auc', gene_set_name=enhan$X, pr=vasari_pr, roc=vasari_roc,
#                               gs_info=gs_info, plot_col_names=labels, top_n=30, type='enhancement',
#                               print=T, h=6.5, w=12, pgratio=c(.8,1), fs=8, both=T, gsname=paste0('enhan'),
#                               clus_row=T, clus_col=T, tags=enhan, show_rows=T, treeheight_row=10,
#                         gs_types=c('canonical','reactome','kegg','biocarta'))

these_gs <- c('canonical','kegg','biocarta','reactome', 'GO', 'motif')

cols <- c(ppalette[-2], mpalette[c(1,4,7)], 'white', ppalette, mpalette)

x <- lapply(diehn_enhan, get_search)
names(x) <- diehn_enhan
x <- bind_rows(x)
x <- x[!(duplicated(x$X)),] # if using as tags, need to remove duplicates
cat('total gs: ', nrow(x))
gene_tag_color <- cols[1:length(diehn_enhan)]; names(gene_tag_color) <- diehn_enhan
diehn_enhan_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                              gs_info=gs_info, plot_col_names=labels, top_n=20, type='enhancing',
                              print=T, h=6, w=11, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('diehn_enhan'),
                              clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, 
                              gs_types=these_gs)

x <- lapply(diehn_infil, get_search); x <- bind_rows(x)
gene_tag_color <- cols[1:length(diehn_infil)]; names(gene_tag_color) <- diehn_infil
diehn_infil_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                            gs_info=gs_info, plot_col_names=labels, top_n=10, type='infiltrative',
                            print=T, h=6, w=12, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('diehn_infil'),
                            clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10)

x <- lapply(gev_necro, get_search); x <- bind_rows(x)
gene_tag_color <- cols[1:length(gev_necro)]; names(gene_tag_color) <- gev_necro
gev_necro_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                            gs_info=gs_info, plot_col_names=labels, top_n=10, type='necrosis',
                            print=T, h=6, w=12, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('gev_necro'),
                            clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10)

x <- lapply(jam_enhan_genes, get_search); x <- bind_rows(x)
x <- rbind(x, data.frame(X=jam_enhan_paths, query='named'))
x <- x[!(duplicated(x$X)),]
gene_tag_color <- cols[1:length(jam_enhan_genes)]; names(gene_tag_color) <- jam_enhan_genes
jam_enhan <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                          gs_info=gs_info, plot_col_names=labels, top_n=20, type='enhancing',
                          print=T, h=6, w=11, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('jam_enhan_all'),
                          clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, gs_types=these_gs)

x <- data.frame(X=jam_enhan_paths, query='named')
gene_tag_color <- cols[1]; names(gene_tag_color) <- 'named'
jam_enhan_1 <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                  gs_info=gs_info, plot_col_names=labels, top_n=10, type='enhancing',
                  print=T, h=5, w=11, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('jam_enhan_biocarta'),
                  clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10,  gs_types=c('biocarta'))

x <- lapply(jam_enhan_genes, get_search); x <- bind_rows(x)
x <- x[!(duplicated(x$X)),]
gene_tag_color <- cols[1:length(jam_enhan_genes)]; names(gene_tag_color) <- jam_enhan_genes
jam_enhan_2 <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                   gs_info=gs_info, plot_col_names=labels, top_n=10, type='enhancing',
                   print=T, h=6, w=11, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('jam_enhan'),
                   clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10,  gs_types=these_gs)


x <- data.frame(X=jam_necro_paths, query='named')
gene_tag_color <- cols[1]; names(gene_tag_color) <- 'named'
jam_necro_1 <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                   gs_info=gs_info, plot_col_names=labels, top_n=10, type='necrosis',
                   print=T, h=8, w=12, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('jam_necro_biocarta'),
                   clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10,  gs_types=c('biocarta'))

x <- lapply(jam_necro_genes, get_search); x <- bind_rows(x)
x <- x[!(duplicated(x$X)),]
gene_tag_color <- cols[1:length(jam_necro_genes)]; names(gene_tag_color) <- jam_necro_genes
jam_necro_2 <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                   gs_info=gs_info, plot_col_names=labels, top_n=10, type='necrosis',
                   print=T, h=8, w=12, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('jam_necro'),
                   clus_row=T, clus_col=F, tags=x, show_rows=T, gs_types=these_gs)

x <- lapply(gut_enhan, get_search)
x <- bind_rows(x)
gene_tag_color <- cols[1:length(gut_enhan)]; names(gene_tag_color) <- gut_enhan
gut_enhan_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                            gs_info=gs_info, plot_col_names=labels, top_n=10, type='enhancing',
                            print=T, h=6, w=11, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('gut_enhan'),
                            clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, gs_types=these_gs)

x <- lapply(gut_necro, get_search); x <- bind_rows(x)
x <- x[!(duplicated(x$X)),]
gene_tag_color <- cols[1:length(gut_necro)]; names(gene_tag_color) <- gut_necro
gut_necro_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                            gs_info=gs_info, plot_col_names=labels, top_n=10, type='necrosis',
                            print=T, h=6, w=11, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('gut_necro'),
                            clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10,  gs_types=these_gs)

x <- lapply(colen_invas, get_search); x <- bind_rows(x)
x <- x[!(duplicated(x$X)),]
gene_tag_color <- cols[1:length(colen_invas)]; names(gene_tag_color) <- colen_invas
colen_invas_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                              gs_info=gs_info, plot_col_names=labels, top_n=10, type='infiltrative',
                              print=T, h=8, w=12, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('colen_invas'),
                              clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10,  gs_types=these_gs)

x <- lapply(zinn_invas_up, get_search); x <- bind_rows(x)
x <- x[!(duplicated(x$X)),]
gene_tag_color <- cols[1:length(zinn_invas_up)]; names(gene_tag_color) <- zinn_invas_up
zinn_edema_up_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                              gs_info=gs_info, plot_col_names=labels, top_n=10, type='edema',
                              print=T, h=8, w=12, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('zinn_edema_up'),
                              clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10,  gs_types=these_gs)

x <- lapply(zinn_invas_dn, get_search); x <- bind_rows(x)
x <- x[!(duplicated(x$X)),]
gene_tag_color <- cols[1:length(zinn_invas_dn)]; names(gene_tag_color) <- zinn_invas_dn
zinn_edema_dn_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                                gs_info=gs_info, plot_col_names=labels, top_n=10, type='edema',
                                print=T, h=8, w=12, pgratio=c(.8,1), fs=9, both=T, gsname=paste0('zinn_edema_dn'),
                                clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10,  gs_types=these_gs)

# GENOMIC gene sets  -----
# gene_tag_color <- RColorBrewer::brewer.pal(name='Set3',n=12)
gene_tag_color <- c(ppalette, mpalette[c(1,4,7)], 'white', ppalette, mpalette)[1:length(gen)]
names(gene_tag_color) <- gen
gene_tag_color[['CDKN2A']] <- 'white'
gene_tag_color[['EGFR']] <- '#5e35b1'

x <- lapply(gen, get_search)
names(x) <- gen
x <- bind_rows(x)
x <- x[!(duplicated(x$X)),] # if using as tags, need to remove duplicates
cat('total gs: ', nrow(x))
genomics_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                              gs_info=gs_info, plot_col_names=labels, top_n=10, type='among_each',
                              print=T, h=6, w=11, pgratio=c(.82,1), fs=10, both=T, gsname=paste0('genomics'),
                              clus_row=T, clus_col=T, tags=x, show_rows=T, treeheight_row=10, 
                              gs_types=these_gs)
dev.print(png, file.path(fdd, paste0('lit_parsons_mclendon','.png')), res=300, height=15, width=12, units="in")
View(vasari_pr[vasari_pr$X %in% x$X & vasari_pr$name %in% c('GO','reactome'),])
View(vasari_roc[vasari_roc$X %in% x$X & vasari_pr$name %in% c('GO','reactome'),])

these_gs <- c('canonical','kegg','biocarta','reactome', 'GO', 'chromosome')
genomics_enhan_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                           gs_info=gs_info, plot_col_names=labels, top_n=15, type='enhancing',
                           print=T, h=6, w=11, pgratio=c(.82,1), fs=9, both=F, gsname=paste0('genomics'),
                           clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, 
                           gs_types=these_gs, title='enhancing AP')
genomics_edema_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                                 gs_info=gs_info, plot_col_names=labels, top_n=15, type='edema',
                                 print=T, h=6, w=11, pgratio=c(.82,1), fs=9, both=F, gsname=paste0('genomics'),
                                 clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, 
                                 gs_types=these_gs, title='edema AP')
genomics_ncet_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                                 gs_info=gs_info, plot_col_names=labels, top_n=15, type='nCET',
                                 print=T, h=6, w=11, pgratio=c(.82,1), fs=9, both=F, gsname=paste0('genomics'),
                                 clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, 
                                 gs_types=these_gs, title='nCET AP')
genomics_focal_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                                 gs_info=gs_info, plot_col_names=labels, top_n=15, type='focal',
                                 print=T, h=6, w=11, pgratio=c(.82,1), fs=9, both=F, gsname=paste0('genomics'),
                                 clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, 
                                 gs_types=these_gs, title='focal AP')
genomics_infil_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                                 gs_info=gs_info, plot_col_names=labels, top_n=15, type='infiltrative',
                                 print=T, h=6, w=11, pgratio=c(.82,1), fs=9, both=F, gsname=paste0('genomics'),
                                 clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, 
                                 gs_types=these_gs, title='infiltrative AP')
genomics_necro_p <- get_heatmaps(metric='pr_auc', gene_set_name=x$X, pr=vasari_pr, roc=vasari_roc, 
                                 gs_info=gs_info, plot_col_names=labels, top_n=15, type='necrosis',
                                 print=T, h=6, w=11, pgratio=c(.82,1), fs=9, both=F, gsname=paste0('genomics'),
                                 clus_row=T, clus_col=F, tags=x, show_rows=T, treeheight_row=10, 
                                 gs_types=these_gs, title='necrosis AP')

# combine gen plots
plot_grid(genomics_enhan_p[[1]][[4]], genomics_edema_p[[1]][[4]],
          genomics_ncet_p[[1]][[4]], genomics_focal_p[[1]][[4]],genomics_infil_p[[1]][[4]], genomics_necro_p[[1]][[4]],
          nrow=3, labels=c('a','b','c','d','e','f',label_size=20, vjust=1))
dev.print(png, file.path(fdd, paste0('lit_parsons_mclendon_individual','.png')), res=300, height=16, width=13, units="in")


## combine enhan plots -----
p1 <- plot_grid(diehn_enhan_p[[1]][[4]], diehn_enhan_p[[2]][[4]], 
                NULL, NULL,  NULL, NULL,nrow=3, rel_widths=c(.8,1), labels=c('a'), label_size=20, vjust=1)
p2a <- plot_grid(jam_enhan_1[[1]][[4]], jam_enhan_1[[2]][[4]], rel_widths=c(.8,1), labels=c('b'), label_size=20, vjust=1)
p2b <- plot_grid(jam_enhan_2[[1]][[4]], jam_enhan_2[[2]][[4]], rel_widths=c(.8,1), labels=c('c'), label_size=20, vjust=1)
p3 <- plot_grid(gut_enhan_p[[1]][[4]], gut_enhan_p[[2]][[4]], rel_widths=c(.8,1), labels=c('d'), label_size=20, vjust=1)

ggdraw() + 
  draw_plot(p1, 0, .05, 1, .95) +
  draw_plot(p2a, 0,.48, 1,.2) +
  draw_plot(p2b, 0,.24, 1,.24) +
  draw_plot(p3, 0,0, 1, .24)
dev.print(png, file.path(fdd, paste0('lit_enhan','.png')), res=300, height=16, width=12, units="in")

p1 <- plot_grid(diehn_enhan_p[[1]][[4]], diehn_enhan_p[[2]][[4]], 
                # NULL, NULL,  
                NULL, NULL,
                nrow=3, rel_widths=c(.8,1), labels=c('a'), label_size=20, vjust=1)
# p2a <- plot_grid(jam_enhan_1[[1]][[4]], jam_enhan_1[[2]][[4]], rel_widths=c(.8,1), labels=c('b'), label_size=20, vjust=1)
p2b <- plot_grid(jam_enhan_2[[1]][[4]], jam_enhan_2[[2]][[4]], rel_widths=c(.8,1), labels=c('c'), label_size=20, vjust=1)
p3 <- plot_grid(gut_enhan_p[[1]][[4]], gut_enhan_p[[2]][[4]], rel_widths=c(.8,1), labels=c('d'), label_size=20, vjust=1)

ggdraw() + 
  draw_plot(p1, 0, -.1, 1, 1.1) +
  draw_plot(p2b, 0,.32, 1,.3) +
  draw_plot(p3, 0,0, 1,.3)
dev.print(png, file.path(fdd, paste0('lit_enhan_smaller','.png')), res=300, height=12, width=13, units="in")

## combine necro plots  -----
p1 <- plot_grid(gev_necro_p[[1]][[4]], gev_necro_p[[2]][[4]], 
                NULL, NULL,  NULL, NULL,nrow=3, rel_widths=c(.8,1), labels=c('a'), label_size=20, vjust=1)
p2a <- plot_grid(jam_necro_1[[1]][[4]], jam_necro_1[[2]][[4]], rel_widths=c(.8,1), labels=c('b'), label_size=20, vjust=1)
p2b <- plot_grid(jam_necro_2[[1]][[4]], jam_necro_2[[2]][[4]], rel_widths=c(.8,1), labels=c('c'), label_size=20, vjust=1)
p3 <- plot_grid(gut_necro_p[[1]][[4]], gut_necro_p[[2]][[4]], rel_widths=c(.8,1), labels=c('d'), label_size=20, vjust=1)

ggdraw() + 
  draw_plot(p1, 0, 0, 1, 1) +
  draw_plot(p2a, 0,.53, 1,.13) +
  draw_plot(p2b, 0,.3, 1,.23) +
  draw_plot(p3, 0,0, 1,.3)
dev.print(png, file.path(fdd, paste0('lit_necro','.png')), res=300, height=15, width=12, units="in")

## combine invasion plots  -----
p1 <- plot_grid(diehn_infil_p[[1]][[4]], diehn_infil_p[[2]][[4]], 
                rel_widths=c(.8,1), labels=c('b'), label_size=20, vjust=1)
p2 <- plot_grid(colen_invas_p[[1]][[4]], colen_invas_p[[2]][[4]], rel_widths=c(.8,1), labels=c('a'), label_size=20, vjust=1)

ggdraw() + 
  draw_plot(p1, 0, 0, 1, .5) +
  draw_plot(p2, 0, .5, 1, .5)
dev.print(png, file.path(fdd, paste0('lit_infil','.png')), res=300, height=10, width=13, units="in")

## combine edema plots  -----
p1 <- plot_grid(zinn_edema_up_p[[1]][[4]], zinn_edema_up_p[[2]][[4]], 
                rel_widths=c(.8,1), labels=c('b'), label_size=20, vjust=1)
p2 <- plot_grid(zinn_edema_dn_p[[1]][[4]], zinn_edema_dn_p[[2]][[4]], rel_widths=c(.8,1), labels=c('a'), label_size=20, vjust=1)

ggdraw() + 
  draw_plot(p1, 0, 0, 1, .5) +
  draw_plot(p2, 0, .5, 1, .5)
dev.print(png, file.path(fdd, paste0('lit_edema','.png')), res=300, height=8, width=12, units="in")


# curve scores vs gene size and coverage  -------

# attch single gene info to gs_info
colnames(gs_info)
sg_info <- data.frame('X'=gids,'cover'=rep(100, length(gids)), 'size'=rep(1, length(gids)))
sg_info$'gene set' <- 'single'

gs_info <- rbind(gs_info, sg_info)
gs_info$size_1k <- gs_info$size/1000


# overall, gene sets combined

vasari <- vasari_pr
vasari <- vasari[, c(labels, 'X')]
vasari <- melt(vasari, variable.name='label')
vasari$score <- 'precision'
r <- vasari_roc
r <- r[, c(labels, 'X')]
r <- melt(r, variable.name='label')
r$score <- 'roc'
vasari <- rbind(vasari, r)
vasari <- left_join(vasari, gs_info[, c('X','cover','size_1k', 'gene set')])
vasari$size_1k <- vasari$size_1k*10

ggplot(na.omit(vasari)) +
  geom_point(aes(x=size_1k, y=value), alpha=0.3, size=1)+ 
  theme_bw() + ylab('performance') + xlab('gene set size (hundreds)') +
  facet_grid(label~score, scales='free_x') 
dev.print(png, file.path(fdd, paste0('gene_set_size_vs_score_all','.png')), res=300, height=10, width=9, units="in")

ggplot(na.omit(vasari)) +
  geom_point(aes(x=cover, y=value), alpha=0.5, size=1)+ 
  theme_bw() + ylab('performance') + xlab('gene set coverage in transcriptome (%)') +
  facet_grid(label~score, scales='free_x') 
dev.print(png, file.path(fdd, paste0('gene_set_coverage_vs_score_all','.png')), res=300, height=10, width=9, units="in")


# for each gene set
vasari <- vasari_pr
vasari <- vasari[, c(labels, 'X')]
vasari <- melt(vasari, variable.name='label', value.name='precision')
vasari <- left_join(vasari, gs_info[, c('X','cover','size_1k', 'gene set')])
vasari$size_1k <- vasari$size_1k*10

ggplot(na.omit(vasari)) +
  geom_point(aes(x=size_1k, y=precision, color=`gene set`), alpha=0.3, show_guide=FALSE, size=1)+ 
  scale_y_continuous(limits=c(.2,1)) +
  theme_bw() + xlab('gene set size (hundreds)') +
  facet_grid(label~`gene set`) 
dev.print(png, file.path(fdd, paste0('gene_set_size_vs_precision','.png')), res=300, height=6, width=12, units="in")

ggplot(na.omit(vasari)) +
  geom_point(aes(x=cover, y=precision, color=`gene set`), alpha=0.3, show_guide=FALSE, size=1)+ 
  theme_bw() + xlab('gene set coverage in transcriptome (%)') +
  facet_grid(label~`gene set`) 
dev.print(png, file.path(fdd, paste0('gene_set_coverage_vs_precision','.png')), res=300, height=6, width=12, units="in")


r <- vasari_roc
r <- r[, c(labels, 'X')]
r <- melt(r, variable.name='label', value.name='roc')
r <- left_join(r, gs_info[, c('X','cover','size_1k', 'gene set')])
r$size_1k <- r$size_1k*10

ggplot(na.omit(r)) +
  geom_point(aes(x=size_1k, y=roc, color=`gene set`), alpha=0.3, show_guide=FALSE, size=1) + 
  scale_y_continuous(limits=c(.3,1)) + 
  theme_bw() + xlab('gene set size (hundreds)') +
  facet_grid(label~`gene set`)   
dev.print(png, file.path(fdd, paste0('gene_set_size_vs_roc','.png')), res=300, height=6, width=12, units="in")


ggplot(na.omit(r)) +
  geom_point(aes(x=cover, y=roc, color=`gene set`), alpha=0.3, show_guide=FALSE, size=1)+ 
  theme_bw() +  xlab('gene set coverage in transcriptome (%)') +
  facet_grid(label~`gene set`)   
dev.print(png, file.path(fdd, paste0('gene_set_coverage_vs_roc','.png')), res=300, height=6, width=12, units="in")

