
library(doParallel)
library(nortest)
library(poLCA)
library(reshape)
library(reshape2)

#####################
# create a formula based on some columns
create_formula <- function(columns) {
    f <- as.formula(paste('cbind(', paste(columns, collapse=', '), ') ~ 1'))
    return(f)
}

#####################
# for the LCA, it is good to fix the variables in case we lose a level while bootstrapping:
fix_categorical <- function(df, columns_to_fix=NULL) {
    if(is.null(columns_to_fix)) {
        columns_to_fix <- names(Filter(is.factor, df))
    }
    for(i in c(columns_to_fix)) { 
        if(i %in% colnames(df)) {
            df[[i]] <- as.factor(as.character(df[[i]]))
            #print(paste(i, paste(levels(df[[i]]), collapse=' ')))
        }
    }
    return(df)
}


#####################
# run lca in a range of groups and seeds (optional)

run_LCA <- function(df, formula, groups=2:7, seeds=1:5, graphs=FALSE, maxiter=5e3, nrep=1) {
  seed <- vector()
  nclasses <- vector()
  bic <- vector()
  aic <- vector()
  llik <- vector()
  models <- list()
  if(is.null(seeds)) {
    s <- 0
    for(g in groups) {
      i <- poLCA(formula,
                df,
                nclass=g,
                maxiter=maxiter,
                tol=1e-10,
                na.rm=FALSE,
                nrep=nrep,
                verbose=TRUE,
                calc.se=TRUE,
                graphs=graphs)
      
      seed <- append(seed, s)
      nclasses <- append(nclasses, g)
      bic <- append(bic, i$bic)
      aic <- append(aic, i$aic)
      llik <- append(llik, i$llik)
      models[[length(seed)]] <- i
    }
  } else {
    for(s in seeds) {
      set.seed(s)
      for(g in groups) {
        i <- poLCA(formula,
                   df,
                   nclass=g,
                   maxiter=maxiter,
                   tol=1e-10,
                   na.rm=FALSE,
                   nrep=nrep,
                   verbose=TRUE,
                   calc.se=TRUE,
                   graphs=graphs)
        
        seed <- append(seed, s)
        nclasses <- append(nclasses, g)
        bic <- append(bic, i$bic)
        aic <- append(aic, i$aic)
        llik <- append(llik, i$llik)
        models[[length(seed)]] <- i
      }
    }
  }
  rdf <- data.frame(seed,
                    nclasses, 
                    bic,
                    aic,
                    llik)
  ret <- list("rdr" = rdf, "models" = models)
  return(ret)
}

#####################
# run lca in a range of groups and seeds (optional)
run_LCA_parallel <- function(df, formula, groups=2:7, seeds=1:5, graphs=FALSE) {
  warning('use operate_LCA and loop_operation_LCA for this analysis')

  seed <- vector()
  nclasses <- vector()
  bic <- vector()
  aic <- vector()
  llik <- vector()
  models <- list()
  if(is.null(seeds)) {
    s <- 0
    ret <- foreach(i=groups, .packages='poLCA') %dopar% poLCA(formula,
                                                              df,
                                                              nclass=i,
                                                              maxiter=5000,
                                                              tol=1e-10,
                                                              na.rm=FALSE,
                                                              nrep=5,
                                                              verbose=TRUE,
                                                              calc.se=TRUE,
                                                              graphs=graphs)
    for(i in ret) {
      seed <- append(seed, s)
      nclasses <- append(nclasses, max(i$predclass))
      bic <- append(bic, i$bic)
      aic <- append(aic, i$aic)
      llik <- append(llik, i$llik)
      models[[length(seed)]] <- i
    }
    
  } else {
    for(s in seeds) {
        set.seed(s)
        ret <- foreach(i=groups, .packages='poLCA') %dopar% poLCA(formula,
                                                                  df,
                                                                  nclass=i,
                                                                  maxiter=5000,
                                                                  tol=1e-10,
                                                                  na.rm=FALSE,
                                                                  nrep=5,
                                                                  verbose=TRUE,
                                                                  calc.se=TRUE,
                                                                  graphs=graphs)
        for(i in ret) {
            seed <- append(seed, s)
            nclasses <- append(nclasses, max(i$predclass))
            bic <- append(bic, i$bic)
            aic <- append(aic, i$aic)
            llik <- append(llik, i$llik)
            models[[length(seed)]] <- i
        }
    }
  }
  rdf <- data.frame(seed,
               nclasses, 
               bic,
               aic,
               llik)
  ret <- list("rdr" = rdf, "models" = models)
  return(ret)
}

#####################
# plot the evolution of the network
plot_performance_lca_clustering <- function(rdf, metric='bic') {
    g <- ggplot(rdf) + 
        geom_boxplot(aes_string(group='nclasses', x='nclasses', y=metric)) + 
        scale_x_discrete(limits=seq(min(rdf$nclasses), max(rdf$nclasses)))    
    return(g)
}

#####################
# run LCA with some settings
run_single_LCA <- function(df, formula, groups, seed, maxiter=1e4, nrep=5) {
    set.seed(seed)
    ret <- poLCA(formula,
                 df,
                 nclass=groups,
                 maxiter=maxiter,
                 tol=1e-10,
                 na.rm=FALSE,
                 nrep=nrep,
                 verbose=TRUE,
                 calc.se=TRUE,
                 graphs=TRUE)
    return(ret)
}

#####################
# plot the polca clusters variables under a stacked bar with the different probabilities

plot_tracks <- function(model) {
  # from: https://statistics.ohlsen-web.de/latent-class-analysis-polca/
  lcmodel <- reshape2::melt(model$probs, level=2)
  zp1 <- ggplot(lcmodel,aes(x = L2, y = value, fill = Var2))
  zp1 <- zp1 + geom_bar(stat = "identity", position = "stack")
  zp1 <- zp1 + facet_grid(Var1 ~ .) 
  zp1 <- zp1 + scale_fill_brewer(type="seq", palette="Greys") + theme_bw()
  zp1 <- zp1 + labs(x = "",y="", fill ="")
  zp1 <- zp1 + theme( axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),                    
                      panel.grid.major.y=element_blank()) + theme(axis.text.x = element_text(angle = 90))
  zp1 <- zp1 + guides(fill = guide_legend(reverse=TRUE))
  return(zp1)
}

#####################
# functions for analysis
operate_LCA <- function(current_seed, groups, df, formula, FILE_FORMAT, nrep=5, redo=FALSE) {
  output_file <- paste0(FILE_FORMAT, max(groups), '_', current_seed, '_.RDS')
  if(file.exists(output_file) & !redo) {
    return(NULL)
  }
  
  set.seed(current_seed)

  ind <- sample.int(nrow(df), size=nrow(df), replace=T)
  fT <- run_LCA(fix_categorical(df[ind, ]), formula, groups=groups, seeds=NULL, nrep=nrep)
  
  saveRDS(fT, output_file)
}


loop_operation_LCA <- function(seeds, groups, df, formula, FILE_FORMAT, nrep=5, redo=FALSE, cores=12) {
  library(parallel)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  foreach(i = seeds, 
          .combine = 'c', 
          .inorder = FALSE, 
          .export = c("operate_LCA", "fix_categorical", "run_LCA"), 
          .packages = c("poLCA", "doParallel")) %dopar% {
            operate_LCA(i, groups, df, formula, FILE_FORMAT, nrep=nrep, redo=redo)
          }
  
  parallel::stopCluster(cl)
}

load_results_LCA <- function(seeds, groups, FILE_FORMAT) {
  final <- NULL
  for(i in seeds) {
    f <- paste0(FILE_FORMAT, max(groups), '_', i,'_.RDS')
    c <- readRDS(f)$rdr
    c$seed <- i
    if(is.null(final)) {
      final <- c
    } else {
      final <- rbind(final, c)
    }
  }
  return(final)
}
