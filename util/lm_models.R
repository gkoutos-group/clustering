library(ModelMetrics)

predict_aucs <- function(model, output_var, predicted_df_train, predicted_df_test) {
    prediction_train <- predict(model, new_data=predicted_df_train)
    prediction_test <- predict(model, new_data=predicted_df_test)

    auc_train <- auc(actual=predicted_df_train[[output_var]], predicted=prediction_train)
    auc_test <- auc(actual=predicted_df_test[[output_var]], predicted=prediction_test)
    
    return(list(train=auc_train, test=auc_test))
}

prepare_OR_table <- function(x) {
    return (cbind(OR=exp(coef(x)), 
                 CI_low=exp(summary(x)$coefficients[, 1] - 1.96*summary(x)$coefficients[, 2]), 
                 CI_high=exp(summary(x)$coefficients[, 1] + 1.96*summary(x)$coefficients[, 2])))
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

lm_models <- function(predicted_df, variables_for_model, case="cluster", col_clusters="predclass", seed=123, train_test_split=0.6) {
    # predicted_df: the dataset with
    # col_clusters: the column indicating the different clusters
    # seed: randomized split of the dataset
    # train_test_split value between 0.01 and 0.99 for splitting the data
    
    f1 <- paste0(variables_for_model, collapse=' + ')
    
    if((train_test_split < 0.01) | (train_test_split > 0.99)) {
        stop("train_test_split must be between 0.01 and 0.99")
    }
    
    set.seed(seed)
    train_ind <- sample(seq_len(nrow(predicted_df)), size=floor(train_test_split * nrow(predicted_df)))

    predicted_df_train <- predicted_df[train_ind,]
    predicted_df_test <- predicted_df[-train_ind,]

    for(i in levels(predicted_df$predclass)) {
        class_is <- paste0('class_', i)
        predicted_df_train[[class_is]] <- as.numeric(predicted_df_train$predclass == i)
        predicted_df_test[[class_is]] <- as.numeric(predicted_df_test$predclass == i)
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
        aucs <- predict_aucs(m1, class_is, predicted_df_train, predicted_df_test)
        or_t1 <- complete_OR_table(m1, class_is, case)
        add_info(case, aucs$train, aucs$test, class_is, or_t1)
    }

    output_results <- cbind(case=output_cases, train=output_auc_train, test=output_auc_test, class=output_class)
    return(list(scores=output_results, or_tables=or_tables))
}
