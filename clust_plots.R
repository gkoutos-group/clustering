

library(ggplot2)


#####################
# simple plot with the number of cases depending on the predicted class
plot_class_distrib <- function(df) {
  g <- ggplot(df, aes(predclass)) + geom_bar()
  return(g)
}

#####################
# plot the values grouped together (stacked bar)
plot_grouped <- function(df, columns_to_test, model, classvar='predclass') {
  classid <- vector()
  comorbidity <- vector()
  probability <- vector()
  #c_to_test <- colnames(df)[grepl('comorbidity.', colnames(df))]
  c_to_test <- columns_to_test
  for(i in c_to_test) {
    comorbidity <- append(comorbidity, i)
    #print(ret$probs[[i]])
    for(j in 1:length(levels(df[[classvar]]))) {
      classid <- append(classid, j)
      probability <- append(probability, model$probs[[i]][j, 2])
    }
  }
  tdf <- data.frame(classid, comorbidity, probability)
  tdf$classid <- as.factor(tdf$classid)
  g <- ggplot(tdf, aes(classid)) + geom_bar(aes(fill = comorbidity, weight = probability))
  return(g)
}

#####################
# for each column plot the distributions
plot_fraction_each <- function(df, columns_to_test, classvar='predclass') {
  #c_to_test <- colnames(df)[grepl('comorbidity.', colnames(df))]
  c_to_test <- columns_to_test
  for(i in c_to_test) {
    tmp_df <- data.frame(ret$probs[[i]])
    colnames(tmp_df) <- c('prob0', 'prob1')
    tmp_df$classes <- as.factor(seq(1, length(levels(df[[classvar]]))))
    #print(tmp_df)
    g <- ggplot(tmp_df, aes(classes)) +
      geom_bar(aes(weight=prob1, fill=prob1)) +
      labs(title = paste("Fraction of comorbidity presence in each class:", i),
           subtitle = "",
           caption = "") +
      ylim(c(0, 1))
    plot(g)
  }    
}

# this function compiles the different useful plots for a complete bootstrapping execution (use with compiled_results)
plot_bootstrapping_results <- function(results) {
    ggplot(results, aes(x=nclasses, group=nclasses, y=bic)) + geom_boxplot() -> boxplts
    results %>% group_by(nclasses) %>% dplyr::summarize(Mean = mean(bic, na.rm=T)) -> overall
    ggplot(results, aes(x=nclasses, group=seed, y=bic, color=seed)) + geom_line() + theme(legend.position = "none") -> bic_lines
    ggplot(results, aes(x=nclasses, y=bic)) + geom_smooth(method='gam') + scale_x_continuous("nclasses", labels=as.character(results$nclasses), breaks=results$nclasses) + theme(legend.position = "none") + theme_minimal() -> bic_smooth

    results %>% group_by(seed) %>% slice(which.min(bic)) -> a
    frequency <- data.frame(table(a$nclasses))
    colnames(frequency) <- c('nclasses', 'frequency')
    ggplot(frequency, aes(x=nclasses, y=frequency)) + geom_bar(stat="identity") -> frequencies
    return(list("boxplts" = boxplts, "overall" = overall, "bic_lines" = bic_lines, "bic_smooth" = bic_smooth, "frequencies" = frequencies))
}
