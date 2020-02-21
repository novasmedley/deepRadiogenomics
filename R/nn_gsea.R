#' Perform GSEA on activations from trained neural net models.
#' Using:
#' a) ranked by correlation (gene to output, need to do)
#' b) ranked by differential expression (? check lit)
#' c) ranked by performance in single gene masking 

library(fgsea)
library(data.table)
library(reactome.db)
library(qvalue)
library(foreach)
library(doParallel)
library(plyr)
library(dplyr)
# dirs
dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular'
fdd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_figs/gsea'

data_dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/with_ae_paper_tabular/nn_retraining_data'
gs_dd <- '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics/gbm_tcga/model_ready'

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

# remove duplicates
gs_info <- split(gs_info, gs_info$`gene set`)
gs_info$canonical <- gs_info$canonical[ !(gs_info$canonical$X %in% c(gs_info$biocarta$X,gs_info$kegg$X, gs_info$reactome$X)), ]
gs_info <- do.call('rbind', gs_info)
rownames(gs_info) <- gs_info$X

#  functions  --------------

gs_names <- c(names(gs_files), 'verhaak','puchalski')
bin_labs <- c('f5','f6','f7','f9','f10','f14')
mul_labs <- c('subtype')
labels <- c(bin_labs, mul_labs)

get_gsea <- function(ranks, rank_df, order_by, gene_set=NA, gmt.file=NA) {
  
  if (is.na(ranks)) {
    rank_df <- rank_df[order(-rank_df[,order_by]), ]
    ranks <- rank_df[, order_by,]
    names(ranks) <- rank_df$X
  } 

  if (!(is.na(gmt.file))) {
    gene_set <- gmtPathways(gmt.file=file.path(gs_dd,'msigdb_v6.2_GMTs', gmt.file))
  } 
  r <- fgsea::fgsea(gene_set, ranks, nperm=10000, maxSize=500, minSize = 15)
  r <- r[order(-r$NES),]
  # r$leadingEdge <- unlist(r$leadingEdge, recursive=F)
  return(r)
}

clean_gene_ids <- function(ids) {
  ids <- gsub("-", ".", ids)
  ids <- gsub("@", ".", ids)
}


rename_classes <- function(c) {
  k <- c
  if (c=='expansive..T1.FLAIR.' | c=='expansive (T1~FLAIR)') k <- 'f10.exp'
  if (c=='infiltrative..T1.much.smaller.than.FLAIR.' | c=='infiltrative (T1 much smaller than FLAIR)') k <- 'f10.mix'
  if (c=='mixed..T1.less.than.FLAIR.' | c=='mixed (T1 less than FLAIR)') k <- 'f10.inf'
  return(k)
}
# correlation --------

get_corr <- function(gs_name, labels) {
  exp_funcs <- c('get_gsea', 'clean_gene_ids', 'mul_labs', 'mapvalues',
                 'rename_classes','bind_rows', 
                 'gs_files', 'gmtPathways', 'gs_dd', 'data_dd', 'verhaak','puchalski')
  
  y_colnames <- c('subtype','f5.proportion.enhancing','f6.proportion.ncet','f7.proportion.necrosis',
                  'f9.multifocal.or.multicentric', 'f10.t1.flair.ratio','f14.proportion.edema')
  names(y_colnames) <- c('subtype','f5','f6','f7','f9','f10','f14')
  
  if (gs_name=='verhaak') gs=verhaak
  if (gs_name=='puchalski') gs=puchalski
  if (gs_name %in% names(gs_files)) gs=gmtPathways(gmt.file=file.path(gs_dd,'msigdb_v6.2_GMTs', gs_files[[gs_name]]))

  parallel::mcaffinity(1:7)
  cl <- makeCluster(7)
  registerDoParallel(cl)
  
  corr <- foreach (i=1:length(labels), .export=exp_funcs ) %dopar% {
    label <- labels[[i]]
    ge <- read.csv(file.path(data_dd, paste0(label, '_original.csv'))) # get genes (scaled differently for each model) and labels
    colnames(ge)[3:12044] <- clean_gene_ids(colnames(ge)[3:12044])  
    
    phenotype <- ge[, y_colnames[label]]
    classes <- levels(phenotype) # less than, and non-focal have class label of 1 during classification
    
    if (label %in% mul_labs) {
      # multiclasses
      pheno <- lapply(seq(1, length(classes)), function(j) {
        indicator <- rep(0, length(classes))
        indicator[j] <- 1
        
        pheno <- mapvalues(phenotype, 
                           from=levels(phenotype),
                           to=indicator) # convert to one-vs-others
        pheno <- as.numeric(as.character(pheno))
      })
      names(pheno) <- classes
    } else {
      # binary
      pheno <- mapvalues(phenotype, 
                         from=classes,
                         to=c(0,1)) # convert to binary ints
      pheno <- as.numeric(as.character(pheno)) 
      pheno <- list(pheno)
      names(pheno) <- label
    }
    
    corrs <- lapply(seq(1,length(pheno)), function(j) {
      ph <- pheno[[j]]
      corr <- apply(ge[, 3:12044], 2, function(k) cor(ph, k))
      corr <- as.data.frame(corr)
      corr$X <- rownames(corr)
      gc <- get_gsea(ranks=NA, rank_df=corr, order_by='corr', gene_set=gs)
      if (nrow(gc) > 0 ) {
        gc$name <- gs_name
        if (label %in% mul_labs) { 
          gc$class <- rename_classes(names(pheno)[j])
        } else {
          gc$class <- label
        }
        gc
      } else {
        NULL
      }
    })
    corrs <- bind_rows(corrs)
  }
  stopCluster(cl)
  return(bind_rows(corr))
}

z <- lapply(gs_names, function(gn) {
  f <- file.path(dd,'gsea', paste('corr_vasari', gn, '.rds', sep='_')); print(f)
  a <- get_corr(gs_name=gn,labels=c(bin_labs))
  if (nrow(a) <1){ print('no results') } else { saveRDS(bind_rows(a), file=f) }
  
  f <- file.path(dd,'gsea', paste('corr_subtype', gn, '.rds', sep='_')); print(f)
  a <- get_corr(gs_name=gn,labels=c('subtype'))
  if (nrow(a) <1){ print('no results') } else { saveRDS(bind_rows(a), file=f) }
})


# single  gsea --------
# class names for reading single masking files
get_class_names <- function(label) {
  if (label=='subtype') {
    cn <- c('Classical','Mesenchymal', 'Neural', 'Proneural')
  }  else {
    cn <- c('')
  }
  print(cn)
  return(cn)
}

get_single <- function(gs_name, labels, metric) {
  parallel::mcaffinity(1:7)
  cl <- makeCluster(7)
  registerDoParallel(cl)
  
  single <- foreach (i=1:length(labels), .export=c('get_gsea', 'clean_gene_ids', 'mul_labs', 'get_class_names',
                                                   'bind_rows',
                                                   'gs_files', 'gmtPathways', 'gs_dd', 'dd', 'verhaak','puchalski') ) %dopar% {
    label <- labels[i]
    if (label %in% mul_labs) {
      cname <- get_class_names(label)
      if (metric=='pr_auc') {
        mname <- '_pr'
      } else {
        mname <- '_roc'
      }
    } else {
      mname <- ''
      if (metric=='pr_auc') {
        cname <- 'pr_auc'
      } else {
        cname <- 'roc_auc'
      }
    }
    if (gs_name=='verhaak') gs=verhaak
    if (gs_name=='puchalski') gs=puchalski
    if (gs_name %in% names(gs_files)) gs=gmtPathways(gmt.file=file.path(gs_dd,'msigdb_v6.2_GMTs', gs_files[[gs_name]]))
      
    s <- read.csv(file.path(dd,'nn_gene_masking', paste0(label, '_single', mname, '_scores.csv')), stringsAsFactors=F) 
    s$X <- clean_gene_ids(s$X)
    
    s_all <- lapply(cname, function(c) {
      s <- get_gsea(ranks=NA, rank_df=s, order_by=c, gene_set=gs)
      s$name <- gs_name
      s$class <- label
      s
    })
    s_all <- bind_rows(s_all)
  }
  stopCluster(cl)
  return(bind_rows(single))
}

z <- lapply(gs_names, function(gn) {
  f <- file.path(dd,'gsea', paste('single_vasari', gn, '.rds', sep='_')); print(f)
  a <- get_single(gs_name=gn,labels=c(bin_labs), metric='pr_auc')
  if (nrow(a) <1){
    print('no results')
  } else {
    saveRDS(bind_rows(a), file=f)
  }
  
  f <- file.path(dd,'gsea', paste('single_subtype', gn, '.rds', sep='_')); print(f)
  a <- get_single(gs_name=gn,labels=c('subtype'), metric='pr_auc')
  if (nrow(a) <1){
    print('no results')
  } else {
    saveRDS(bind_rows(a), file=f)
  }
})




# activations --------
# each file has ONE class activations

get_amax <- function(gs_name, labels) {
  parallel::mcaffinity(1:7)
  cl <- makeCluster(7)
  registerDoParallel(cl)

  single <- foreach (i=1:length(labels), .export=c('get_gsea', 'clean_gene_ids', 'mul_labs',
                                                   'gs_files', 'gmtPathways', 'gs_dd', 'dd',
                                                   'verhaak','puchalski') ) %dopar% {
   label <- labels[i]
   if (label=='subtype') {
     nodes <- as.character(seq(0,3,1))
   } else if (labels=='f10') {
     nodes <- as.character(seq(0,2,1))
   } else {
     nodes <- c('0')
   }
   if (gs_name=='verhaak') gs=verhaak
   if (gs_name=='puchalski') gs=puchalski
   if (gs_name %in% names(gs_files)) gs=gmtPathways(gmt.file=file.path(gs_dd,'msigdb_v6.2_GMTs', gs_files[[gs_name]]))

   a <- lapply(nodes, function(node) {
     amax <- read.csv(file.path(dd,'nn_gene_importance', paste0('activations_', label,'_node_', as.character(node),'.csv')), stringsAsFactors=F)
     amax$X <- clean_gene_ids(amax$X)
     amax <- get_gsea(rank=NA, rank_df=amax, order_by='activations', gene_set=gs)
     amax$name <- gs_name
     amax
   })
 }
  stopCluster(cl)
  return(a)
}

a <- get_amax(gs_name='verhaak',labels='f9')

z <- lapply(gs_names, function(gn) {
  f <- file.path(dd,'gsea', paste('actmax_vasari', gn, '.rds', sep='_')); print(f)
  a <- get_amax(gs_name=gn,labels=c(bin_labs, 'f10'))
  saveRDS(bind_rows(a), file=f)

  f <- file.path(dd,'gsea', paste('actmax_subtype', gn, '.rds', sep='_')); print(f)
  a <- get_amax(gs_name=gn,labels=c('subtype'))
  saveRDS(bind_rows(a), file=f)
})

# saliency --------

# in saliency, each file has ONE class activations
# class names for reading saliency files
get_class_names <- function(label) {
  if (label=='subtype') {
    cn <- c('Classical','Mesenchymal', 'Neural', 'Proneural')
  }  else {
    cn <- c('')
  }
  print(cn)
  return(cn)
}

get_saliency <- function(label, gs, gs_name){
  f <- file.path(dd,'gsea', paste('saliency', label, gs_name, '.rds', sep='_')); print(f)
  
  parallel::mcaffinity(1:7)
  cl <- makeCluster(7)
  registerDoParallel(cl)
  
  cn <- get_class_names(label)
  zp <- 0
  sal_all <- lapply(cn, function(this_class) {
    print(this_class)
    sal_fn <- file.path(dd,'nn_gene_importance', paste0('saliency_patientwise_', label, '_', this_class, '.csv'));print(sal_fn) # saliency data
    sal <- read.csv(sal_fn, stringsAsFactors=F)
    class_labs <- sal[sal$X=='y_label',][,-1]
    
    sal$X <- clean_gene_ids(sal$X)  #prep 
    sal <- sal[sal$X!='y_label',]
    patients <- colnames(sal)[2:ncol(sal)]
    sal[, patients] <- apply(sal[, patients], 2, function(i) as.numeric(i))
    
    # if (label %in% c('subtype')) {
    #   c <- get_real_class_name(this_class); print(c %in% class_labs) # get just the patients with the label
    #   patients <- names(class_labs)[class_labs==c]; length(these)
    # } 
    
    z <- apply(sal[, patients], 2, sum) # check zero saliency
    zp <<- patients[which(z==0)]
    cat('zero sal. patients: ', zp, '\n')
    if (length(zp)!=0) patients <- patients[which(z!=0)] # remove zp patients from list

    start_time <- Sys.time()
    s <- foreach (i=1:length(patients), .export=c('get_gsea', 'gs') ) %dopar% {
      j <- patients[i]
      a <- get_gsea(ranks=NA, rank_df=sal[, c('X',j)], order_by=j, gene_set=gs)[, c('pathway','NES','padj')]
      a$name <- j
      a$class <- this_class
      a
    }
    end_time <- Sys.time(); t2 <- end_time - start_time ; print(t2)
    s <- bind_rows(s)

  })
  s <- list(sal=bind_rows(sal_all), zero_sal=zp)
  saveRDS(s, f)
  
  stopCluster(cl)
  return(s)
}

get_saliency_all <- function(labels, these_gs) {
  lapply(labels, function(label) {
    lapply(these_gs, function(g) {
      print(label)
      print(g)
      if (g=='verhaak') z <- get_saliency(label=label, gs=verhaak, gs_name=g)
      if (g=='puchalski') z <- get_saliency(label=label, gs=puchalski, gs_name=g)
      if (g %in% names(gs_files)) z <- get_saliency(label=label,
                                                    gs=gmtPathways(gmt.file=file.path(gs_dd,'msigdb_v6.2_GMTs', gs_files[[g]])),
                                                    gs_name=g)
      NULL
    })
    NULL
  })
}

names(gs_files) # possible gene sets


these <- c('verhaak','puchalski','hallmark','chromosome','reactome','biocarta','kegg','canonical', 'motif','GO','onco')
these <- c('reactome','biocarta','kegg','canonical', 'motif','GO','onco')
get_saliency_all(labels=c('f5'), these)

these <- c('chem', 'computational','immuno')
get_saliency_all(labels=c(bin_labs, mul_labs), these)
