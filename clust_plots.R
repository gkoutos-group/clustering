

library(ggplot2)


#####################
# simple plot with the number of cases depending on the predicted class
plot_class_distrib <- function(df) {
  g <- ggplot(df, aes(predclass)) + geom_bar()
  return(g)
}

#####################
# plot the values grouped together (stacked bar)
plot_grouped <- function(df, columns_to_test) {
  classid <- vector()
  comorbidity <- vector()
  probability <- vector()
  #c_to_test <- colnames(df)[grepl('comorbidity.', colnames(df))]
  c_to_test <- columns_to_test
  for(i in c_to_test) {
    comorbidity <- append(comorbidity, i)
    #print(ret$probs[[i]])
    for(j in 1:4) {
      classid <- append(classid, j)
      probability <- append(probability, ret$probs[[i]][j, 2])
    }
  }
  tdf <- data.frame(classid, comorbidity, probability)
  tdf$classid <- as.factor(tdf$classid)
  g <- ggplot(tdf, aes(classid)) + geom_bar(aes(fill = comorbidity, weight = probability))
  return(g)
}

#####################
# for each column plot the distributions
plot_fraction_each <- function(df, columns_to_test) {
  #c_to_test <- colnames(df)[grepl('comorbidity.', colnames(df))]
  c_to_test <- columns_to_test
  for(i in c_to_test) {
    tmp_df <- data.frame(ret$probs[[i]])
    colnames(tmp_df) <- c('prob0', 'prob1')
    tmp_df$classes <- as.factor(seq(1, 4))
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
