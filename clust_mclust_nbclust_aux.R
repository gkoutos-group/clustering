
library(mclust)
library(NbClust)

#####################
# execute mclust and save
run_mclust <- function(df, groups, modelNames=NULL, save_to=NULL) {
  mbic <- Mclust(df, G=groups, modelNames=modelNames, verbose=F)
  if(!is.null(save_to)) {
    saveRDS(mbic, file=save_to)
  }
  return(mbic)
}

#####################
# execute nbclust and save
run_nbclust <- function(df, groups, method='kmeans', save_to=NULL) {
  min.nc <- max(min(groups), 2)
  max.nc <- max(groups, min.nc + 2)
  nret <- NbClust(df, method='kmeans', min.nc=min.nc, max.nc=max.nc)
  if(!is.null(save_to)) {
    saveRDS(nret, file=save_to)
  }
  return(nret)
}

#####################
# execute both mclust and nbclust
run_mclust_nbclust <- function(df, groups, save_nbclust=NULL, save_mclust=NULL, save_to=NULL) {
  ret <- list(nbclust=run_nbclust(df, groups, save_to=save_nbclust),
       mclust=run_mclust(df, groups, save_to=save_mclust))
  if(!is.null(save_to)) {
    saveRDS(ret, save_to)
  }
  return(ret)
}

#####################
# operates mclust on a range of different clusters and obtain other metrics not available directly

# from: https://stackoverflow.com/questions/28342653/model-selection-with-aic
get_aic_llik_scores <- function(df, columns_to_test, modelNames='VVI', range_to_test=1:25) {
  index <- vector()
  aics <- vector()
  llik <- vector()
  for(i in range_to_test) {
    index <- append(index, i)
    IC <- Mclust(data=df[, columns_to_test], modelNames=modelNames, G=i)
    aic <- 2*IC$df - 2*IC$loglik
    aics <- append(aics, aic)
    llik <- append(llik, IC$loglik)
    print(paste(i, aic))
  }
  plot(index, aics)
  return(data.frame(index, aics, llik))
}

#####################
# different functions to run multiple tests
operate_MCLUST_NBCLUST <- function(current_seed, groups, df, FILE_FORMAT, redo=FALSE, do_mclust=T, do_nbclust=T) {
  output_file <- paste0(FILE_FORMAT, max(groups), '_', current_seed, '_.RDS')
  if(file.exists(output_file) & !redo) {
    return(NULL)
  }
  
  set.seed(current_seed)
  
  ind <- sample.int(nrow(df), size=nrow(df), replace=T)
  fT <- NULL
  if(do_mclust & do_nbclust) {
    fT <- run_mclust_nbclust(df[ind, ], groups=groups)
  } else if(do_mclust) {
    fT <- run_mclust(df[ind, ], groups=groups)    
  } else if(do_nbclust) {
    fT <- run_nbclust(df[ind, ], groups=groups)
  } else {
    stop('Neither nbclust or mclust selected')
  }
  
  cat(paste("Saved to:", output_file))
  saveRDS(fT, output_file)
}


loop_operation_MCLUST_NBCLUST <- function(seeds, groups, df, FILE_FORMAT, redo=FALSE, do_mclust=T, do_nbclust=T) {
  library(parallel)
  cl <- parallel::makeCluster(12)
  doParallel::registerDoParallel(cl)
  
  
  foreach(i = seeds, 
          .combine = 'c', 
          .inorder = FALSE, 
          .export = c("operate_MCLUST_NBCLUST", "run_mclust_nbclust", "run_mclust", "run_nbclust"), 
          .packages = c("mclust", "NbClust", "doParallel")) %dopar% {
            operate_MCLUST_NBCLUST(i, groups, df, FILE_FORMAT, redo=redo, do_mclust=do_mclust, do_nbclust=do_nbclust)
          }
  
  parallel::stopCluster(cl)
}

load_results_MCLUST_NBCLUST <- function(seeds, groups, FILE_FORMAT, do_mclust=T, do_nbclust=T) {
  # there is a bit of copy and paste in this function
  final_a <- NULL # will store the final result
  final_b <- NULL
  
  result_mclust <- function(a, seed) {
    a$groups <- rownames(a)
    rownames(a) <- NULL
    a$seed <- seed
    if(is.null(final_a)) {
      final_a <- a
    } else {
      final_a <- rbind(final_a, a)
    }
    return(final_a)
  }
  
  result_nbclust <- function(b, seed) {
    b$Metric <- rownames(b)
    rownames(b) <- NULL
    b$seed <- seed
    if(is.null(final_b)) {
      final_b <- b
    } else {
      final_b <- rbind(final_b, b)
    }
    return(final_b)
  }
  
  for(i in seeds) { # iterate through the seeds
    f <- paste0(FILE_FORMAT, max(groups), '_', i,'_.RDS')
    c <- readRDS(f)
    
    if(do_mclust & do_nbclust) { #both cases
      final_a <- result_mclust(data.frame(c$mclust$BIC[,]), i)
      final_b <- result_nbclust(data.frame(t(c$nbclust$Best.nc)), i)
    } else if(do_mclust) { #fix the mclust set
      final_a <- result_mclust(data.frame(c$BIC[,]), i)
    } else if(do_nbclust) { #fix the nbclust set
      final_b <- result_nbclust(data.frame(t(c$Best.nc)), i)
    } else {
      stop('Neither mclust or nbclust selected')
    }
  }
  final <- NULL
  if(do_mclust & do_nbclust) {
    final <- list(mclust=final_a, nbclust=final_b)
  } else if(do_mclust) {
    final <- final_a
  } else if(do_nbclust) {
    final <- final_b
  } else {
    stop('Neither mclust or nbclust selected')
  }
  return(final)
}

