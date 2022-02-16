library(ModelMetrics)
library(pROC)
library(checkmate)
library(MASS)

predict_aucs <- function(model, output_var, predicted_df_train, predicted_df_test) {
    prediction_train <- predict(model, newdata=predicted_df_train, type='response')
    prediction_test <- predict(model, newdata=predicted_df_test, type='response')
    
    auc_train <- auc(predicted_df_train[[output_var]], prediction_train, levels=c(0, 1), direction='<')
    auc_test <- auc(predicted_df_test[[output_var]], prediction_test, levels=c(0, 1), direction='<')
    
    return(list(train=auc_train, test=auc_test))
}

predict_rmse <- function(model, output_var, predicted_df_train, predicted_df_test) {
    prediction_train <- predict(model, newdata=predicted_df_train, type='response')
    prediction_test <- predict(model, newdata=predicted_df_test, type='response')
    
    rmse_train <- rmse(predicted_df_train[[output_var]], prediction_train)
    rmse_test <- rmse(predicted_df_test[[output_var]], prediction_test)
    
    return(list(train=rmse_train, test=rmse_test))
}

prepare_OR_table <- function(x) {
    return (cbind(OR=exp(coef(x)), 
                 CI_low=exp(summary(x)$coefficients[, 1] - 1.96*summary(x)$coefficients[, 2]), 
                 CI_high=exp(summary(x)$coefficients[, 1] + 1.96*summary(x)$coefficients[, 2]),
                 pval=summary(x)$coefficients[, 4]))
}

hr_coefs <- function(df, predictors, tvar, svar) {
    tformula <- as.formula(paste0('Surv(time=', tvar, ', event=', svar, ') ~ ', paste(predictors, collapse=' + ')))
    res.cox <- coxph(tformula, data=df)
    
    scx <- as.data.frame(summary(res.cox)$coefficients)
    scx$HRlow <- exp(scx[, 1] - 1.96*scx[, 3])
    scx$HRhigh <- exp(scx[, 1] + 1.96*scx[, 3])
    return(list(model=res.cox, coeffs=scx))
}

complete_OR_table <- function(x, class_is, case) {
    or_t <- data.frame(prepare_OR_table(x))
    or_t$variable <- rownames(or_t)
    or_t$class <- class_is
    or_t$case <- case
    rownames(or_t) <- NULL
    or_t <- or_t[c('variable', setdiff(colnames(or_t), c('variable')))]
    return(or_t)
}


lm_model_continuous <- function(predicted_df, 
                                variables_for_model, 
                                case="cluster", 
                                output="predclass",
                                seed=123, 
                                train_test_split=0.6,
                                direction=NULL) {
    assert(checkChoice(direction, c('backward', 'forward', 'both'), null.ok=T))
    
    f1 <- paste0(variables_for_model, collapse=' + ')
    
    set.seed(seed)
    train_ind <- sample(seq_len(nrow(predicted_df)), size=floor(train_test_split * nrow(predicted_df)))
    
    predicted_df_train <- predicted_df[train_ind,]
    predicted_df_test <- predicted_df[-train_ind,]

    m1 <- glm(formula(paste0(output, ' ~ ', f1)), family='poisson', data=predicted_df_train)
    if(!is.null(direction)) {
        print('Backward selection process...')
        sdir = stepAIC(m1,
                       direction=direction,
                       trace=T)
        m1 <- sdir
    }
    
    rmses <- predict_rmse(m1, output, predicted_df_train, predicted_df_test)
    return(list(scores=as.data.frame(list(case=case, train=rmses$train, test=rmses$test, output=output)), or_tables=summary(m1)$coefficients))
}

lm_models <- function(predicted_df, 
                      variables_for_model, 
                      case="cluster", 
                      col_clusters="predclass",
                      seed=123, 
                      train_test_split=0.6,
                      direction=NULL) {
    # predicted_df: the dataset with
    # col_clusters: the column indicating the different clusters
    # seed: randomized split of the dataset
    # train_test_split value between 0.01 and 0.99 for splitting the data
    # direction: forward, backward or both - feature selection
    
    assert(checkChoice(direction, c('backward', 'forward', 'both'), null.ok=T))
    
    f1 <- paste0(variables_for_model, collapse=' + ')
    
    if((train_test_split < 0.01) | (train_test_split > 0.99)) {
        stop("train_test_split must be between 0.01 and 0.99")
    }
    
    set.seed(seed)
    train_ind <- sample(seq_len(nrow(predicted_df)), size=floor(train_test_split * nrow(predicted_df)))

    predicted_df_train <- predicted_df[train_ind,]
    predicted_df_test <- predicted_df[-train_ind,]

    for(i in levels(predicted_df[[col_clusters]])) {
        class_is <- paste0('class_', i)
        predicted_df_train[[class_is]] <- as.numeric(predicted_df_train[[col_clusters]] == i)
        predicted_df_test[[class_is]] <- as.numeric(predicted_df_test[[col_clusters]] == i)
    }
    
    or_tables <- NULL
    output_cases <- vector()
    output_auc_train <- vector()
    output_auc_test <- vector()
    output_class <- vector()

    add_info <- function(case, train_auc, test_auc, class_is, or_t) {
        output_cases <<- append(output_cases, case)
        output_auc_train <<- append(output_auc_train, train_auc)
        output_auc_test <<- append(output_auc_test, test_auc)
        output_class <<- append(output_class, class_is)

        if(is.null(or_tables)) {
            or_tables <<- or_t
        } else {
            or_tables <<- rbind(or_tables, or_t)
        }
    }

    for(i in levels(predicted_df[[col_clusters]])) {
        class_is <- paste0('class_', i)

        m1 <- glm(formula(paste0(class_is, ' ~ ', f1)), family=binomial(link='logit'), data=predicted_df_train)
        if(!is.null(direction)) {
            print('Backward selection process...')
            sdir = stepAIC(m1,
                        direction=direction,
                        trace=T)
            m1 <- sdir
        }
        aucs <- predict_aucs(m1, class_is, predicted_df_train, predicted_df_test)
        or_t1 <- complete_OR_table(m1, class_is, case)
        add_info(case, aucs$train, aucs$test, class_is, or_t1)
    }

    output_results <- cbind(case=output_cases, train=output_auc_train, test=output_auc_test, class=output_class)
    return(list(scores=output_results, or_tables=or_tables))
}

# run_for function creates a model and plot auc and ci
# model will be created towards @param outcome
# model will be created on the `train` df in the `dfs` list
# model will be weighted by distribution of values
# @param outcome: outcome variable to be tested against; must be 1 or 0
# @param variables: variables used to created model
# @param predict_on: dataset in `dfs` list used to predict the result
# @param dfs: list of dataframes, `train` and the one pointed by `predict_on` are required
# @return: the predictions made on the new data
lm_balanced <- function(outcome, variables, predict_on, dfs) {
    f <- as.formula(paste0(paste0(outcome, ' ~'), 
                      paste(variables, collapse=' + ', sep=' + ')))
    wt <- 1 + (as.numeric(as.character(dfs[['train']][[outcome]])) * (as.numeric(nrow(dfs[['train']]) / table(dfs[['train']][[outcome]])['1'] - 1)))
    print(f)
    model <- glm(f, 
                 weights = wt,
                 data=dfs[['train']], 
                 family=binomial("logit"))

    tdf <- dfs[[predict_on]]
    pred <- predict(model, newdata=tdf, type='response')

    rocobj_c <- plot.roc(tdf[[outcome]], pred, direction='<', levels=c(0, 1))
    ci.sp.obj_c <- ci.sp(rocobj_c, sensitivities=seq(0, 1, .01), boot.n=100)
    #plot(rocobj_c)
    #plot(ci.sp.obj_nnc, type="shape", col=alpha("green", 0.6), add=T)

    print(rocobj_c$auc)
    print(ci.auc(rocobj_c, method='delong'))
    return(list(model=model, pred=pred))
}
