
library(doParallel)
library(nortest)
library(poLCA)
library(reshape)
library(reshape2)
library(hash)
library(fossil)

#####################
# create a formula based on some columns
create_formula <- function(columns) {
        f <- as.formula(paste('cbind(', paste(columns, collapse=', '), ') ~ 1'))
        return(f)
}

#####################
# for the LCA, it is good to fix the variables in case we lost a level while bootstrapping:
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
    # these vectors will store the whole set of results
    seed <- vector()
    nclasses <- vector()
    bic <- vector()
    aic <- vector()
    llik <- vector()
    models <- list()
        
    if(is.null(seeds)) { # only a single execution
        seeds <- c(0)
    }

    for(s in seeds) { # loop through the seeds
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
    
    # combine the vectors
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
    warning('this function is outdated and should be avoided')

    # these vectors will store the whole set of results
    seed <- vector()
    nclasses <- vector()
    bic <- vector()
    aic <- vector()
    llik <- vector()
    models <- list()
    
    if(is.null(seeds)) { # only a single execution
        seeds <- c(0)
    }
        
    for(s in seeds) { # loop through the seeds
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
            for(i in ret) { # combine the foreach return
                    seed <- append(seed, s)
                    nclasses <- append(nclasses, max(i$predclass))
                    bic <- append(bic, i$bic)
                    aic <- append(aic, i$aic)
                    llik <- append(llik, i$llik)
                    models[[length(seed)]] <- i
            }
    }
    
    # combine the vectors
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
run_single_LCA <- function(df, formula, groups, seed, maxiter=1e4, nrep=5, graphs=TRUE) {
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
                     graphs=graphs)
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
operate_LCA <- function(current_seed, groups, df, formula, FILE_FORMAT, nrep=5, redo=FALSE, sampling_rate=1) {
    output_file <- paste0(FILE_FORMAT, max(groups), '_', current_seed, '_.RDS')
    if(file.exists(output_file) & !redo) {
        return(NULL)
    }
    
    set.seed(current_seed)

    ind <- sample.int(nrow(df), size=as.integer(sampling_rate*nrow(df)), replace=T)
    fT <- run_LCA(fix_categorical(df[ind, ]), formula, groups=groups, seeds=NULL, nrep=nrep)
    
    saveRDS(fT, output_file)
}


loop_operation_LCA <- function(seeds, groups, df, formula, FILE_FORMAT, nrep=5, redo=FALSE, cores=12, rslurm=FALSE, sampling_rate=1, rslurm_qos="castles", rslurm_account="gkoutosg-variant-prediction", rslurm_jobname='clustering') {
    if(rslurm == TRUE) {
        library(rslurm)
        rslurm_operate <- function(i) {
            operate_LCA(i, groups, df, formula, FILE_FORMAT, nrep=nrep, redo=redo, sampling_rate=sampling_rate)
        }
        sjob <- slurm_apply(rslurm_operate,
                            params=data.frame(i=seeds), 
                            add_objects=c('groups', 'df', 'formula', 'FILE_FORMAT', 'nrep', 'redo', 'operate_LCA', 'fix_categorical', 'run_LCA', 'create_formula', 'var_comorb', 'sampling_rate'), 
                            nodes=cores,
                            cpus_per_node=1, 
                            pkgs=c('poLCA'),
                            slurm_options = list(time = "8:00:00",
                                                 qos =rslurm_qos, account=rslurm_account), 
                            jobname=rslurm_jobname)
        cleanup_files(sjob)
    } else {
        library(parallel)
        cl <- parallel::makeCluster(cores, outfile='temp_op_LCA_parallel.out')
        doParallel::registerDoParallel(cl)
        
        foreach(i = seeds, 
                        .combine = 'c', 
                        .inorder = FALSE, 
                        .export = c("operate_LCA", "fix_categorical", "run_LCA"), 
                        .packages = c("poLCA", "doParallel")) %dopar% {
                            operate_LCA(i, groups, df, formula, FILE_FORMAT, nrep=nrep, redo=redo, sampling_rate=sampling_rate)
                        }
        
        parallel::stopCluster(cl)
    }
}

rand_index_pairwise <- function(models) {
    modelA <- models[[1]]
    modelB <- models[[2]]
    
    # prepare the data from models
    mA <- unique(data.frame(ids=rownames(modelA$x), predclass=modelA$predclass))
    mB <- unique(data.frame(ids=rownames(modelB$x), predclass=modelB$predclass))
    
    # obtain the rows in common
    rows_to_use <- intersect(mA$ids, mB$ids)
    mA <- mA[mA$ids %in% rows_to_use,]
    mB <- mB[mB$ids %in% rows_to_use,]
    
    # just to be sure the comparison is against the same ids
    if((nrow(mA) == 0) | (nrow(mB) == 0)) {
        return(0)
    }
    mA <- mA[with(mA, order(ids)),]
    mB <- mB[with(mB, order(ids)),]
        
    # calculate and return index
    ri <- rand.index(mA$predclass, mB$predclass)
    return(ri)
}

# gets the rand index for the results
compile_rindex_results <- function(best_models, parallel_cores) {
    compiled_results <- list()
    # go through the list of best models with N classes
    for(i in keys(best_models)) {
        case <- values(best_models, i)
        if(parallel_cores == 1){
            # compare each pair collecting the rand-index score
            compiled_results[[i]] <- apply(combn(2:length(case), 2), 2, function(x) rand_index_pairwise(case[x]))
        } else {
            if(.Platform$OS.type == "unix") {
                library(parallel)
            } else {
                library(parallelsugar)
            }
            #cl <- parallel::makeCluster(detectCores()-1)
            #doParallel::registerDoParallel(cl)

            combinations <- t(combn(2:length(case), 2))
            ri_values <- mclapply(1:nrow(combinations), function(x) rand_index_pairwise(case[combinations[x, ]]), mc.cores=parallel_cores)
            
            ri_values[sapply(ri_values, is.null)] <- NA
                            
            compiled_results[[i]] <- data.frame(first=unlist(combinations[, 1]), second=unlist(combinations[, 2]), value=unlist(ri_values))

            #parallel::stopCluster(cl)
        }
    }
    return(compiled_results)
}

# this function loads the results from bootstrapping
load_results_LCA <- function(seeds, groups, FILE_FORMAT, obtain_rindex=FALSE, default_length=10, optimization_parameter='bic', parallel_cores=6) {
    final <- NULL
    rindex <- NULL
    best_models <- hash()
        
    for(i in seeds) {
        f <- paste0(FILE_FORMAT, max(groups), '_', i,'_.RDS')
        f_content <- readRDS(f)
        c <- f_content$rdr
        c$seed <- i
        if(is.null(final)) {
            final <- c
        } else {
            final <- rbind(final, c)
        }
        
        if(! is.null(optimization_parameter)) {
            best_model <- c[which.min(c[[optimization_parameter]]),]
            this_class <- best_model$nclasses
            this_class_char <- as.character(this_class)
            model <- f_content$models[[this_class - min(c$nclasses) + 1]] # the list might not start with 1
                        
            if(! has.key(this_class_char, best_models)) {
                best_models[this_class_char] <- list()
            }
            stored_models <- values(best_models, this_class_char)
            stored_models[[length(stored_models) + 1]] <- model
            best_models[this_class_char] <- stored_models
        }
    }
    
    if(obtain_rindex) {
        rindex <- compile_rindex_results(best_models,
                                         parallel_cores=parallel_cores)
    } else {
        rindex <- NULL
    }
    
    return(list("final" = final, "rindex" = rindex))
}

warning('make sure that _clustering/clust_info_table.R_ is loaded as well')
# for a number of clusters, execute and save everything
run_and_save <- function(df,
                         formula, 
                         N_clusters, 
                         output_format, # to save the files
                         var_numerical,
                         var_categorical,
                         var_comorb,
                         seed=1, 
                         nrep=1,
                         verbose=TRUE,
                         transform_cat=F,
                         graphs=FALSE) {
    if(verbose) {
        cat('creating model... ')
    }
    model <- run_single_LCA(df, formula, N_clusters, seed=seed, nrep=nrep, graphs=graphs)
    
    predicted_df <- df
    predicted_df$predclass <- model$predclass
    
    if(verbose) {
        cat('preparing tables... ')
    }
    
    if(transform_cat) {
        for(i in var_categorical) {
            predicted_df[[i]] <- as.factor(predicted_df[[i]])
        }
    }
    
    ret <- compile_results_to_xlsx(predicted_df, 
                        continuous_variables=var_numerical, 
                        categorical_variables=var_categorical, 
                        comorbidity_variables=var_comorb,
                        output_file=paste0(output_format, ".xlsx"),
                        subgroup_cases=c(1, 2, 3, 4, 5, 6, 7, 8),
                        positive_class="1",
                        shapiro_threshold=0.05,
                        cname='comorbidities', 
                        cvalue="1",
                        classvar='predclass')
    
    return(predicted_df)
}
