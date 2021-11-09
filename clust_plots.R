

library(ggplot2)
library(plyr)


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
    ggplot(results, aes(x=nclasses, y=bic)) + geom_smooth(method='gam', formula = y ~ s(x, bs='cs', k=5)) + scale_x_continuous("nclasses", labels=as.character(results$nclasses), breaks=results$nclasses) + theme(legend.position = "none") + theme_minimal() -> bic_smooth

    results %>% group_by(seed) %>% slice(which.min(bic)) -> a
    frequency <- data.frame(table(a$nclasses))
    colnames(frequency) <- c('nclasses', 'frequency')
    ggplot(frequency, aes(x=nclasses, y=frequency)) + geom_bar(stat="identity") -> frequencies
    return(list("boxplts" = boxplts, "overall" = overall, "bic_lines" = bic_lines, "bic_smooth" = bic_smooth, "frequencies" = frequencies))
}

plot_heatmap_relative <- function(df, 
                                  groupvar='predclass', 
                                  variables=NULL,
                                  variables_renames=NULL, # in format c(old="New", other_old="New other")
                                  positive_values=c(1, "1", 'Y', 'y', "Yes", "yes"),
                                  palette=NULL, # if null use the colors on palette colors from brewer.pal
                                  palette_colors='RdYlBu',
                                  output_overall_heatmap_file=NULL,
                                  output_r_heatmap_file=NULL,
                                  output_csv_values=NULL
                                  ) {
  # if no variable is specific get them
  if(is.null(variables)) {
    variables <- setdiff(colnames(df), c(groupvar))
    variables <- variables[order(variables)]
  }
  
  # get the color definition
  if(is.null(palette)) {
    palette <- rev(brewer.pal(n=10, name=palette_colors))
  }
  
  ## prepare the main dataset
  ready_for_heat <- df[, c(variables, groupvar)] %>% 
    melt(id=c(groupvar)) %>% 
    group_by(predclass, variable, value) %>% 
    count() %>% 
    group_by(predclass, variable) %>% 
    mutate(percent = n/sum(n))
  
  # reorder elements
  ready_for_heat$variable <- factor(ready_for_heat$variable,
                                    levels=rev(variables))
  
  # rename variables
  if(!is.null(variables_renames)) {
    ready_for_heat$variable <- recode(ready_for_heat$variable, !!!variables_renames)
  }
  
  # only use positive values, first checking them
  ready_for_heat_check <- ready_for_heat[ready_for_heat$value %in% positive_values, ]
  
  # check if all the variables are kept
  if(length(unique(ready_for_heat$variable)) != length(unique(ready_for_heat_check$variable))) {
    stop(paste('Some variables are not going through on the analysis, check if the variables "', 
               paste(setdiff(unique(ready_for_heat$variable), unique(ready_for_heat_check$variable)), sep='", "'),
               '" have their corresponding positive value in "positive_values".'))
  }
  
  # check if there is a single factor option in the list
  if(length(unique(ready_for_heat_check$variable)) != length(ready_for_heat_check$variable)) {
    stop(paste('Some variables are showing more than once, check if the variables "',
               paste(ready_for_heat_check$variable[duplicated(ready_for_heat_check$variable)], sep='", "'),
               '" have a single factor in "positive_values".'))
  }
  
  # checks complete
  ready_for_heat <- ready_for_heat_check
  
  # overall % heatmap
  main_heatmap <- ready_for_heat %>% 
    ggplot(aes_string(x=predclass, y='variable')) + 
    geom_tile(aes(fill=percent)) + 
    theme_minimal() + 
    scale_colour_stepsn(colors = palette, aesthetics = "fill",
                        limits=c(0, 1),
                        breaks=seq(0.1, 0.9, 0.1)) +
    ggtitle(' ') +
    theme(legend.position='none',
          #axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          #axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  ### overall pop
  ready_for_overall <- df[, c(variables)] %>% 
    melt(id=c()) %>% 
    group_by(variable, value) %>% 
    count() %>% 
    group_by(variable) %>% 
    mutate(percent = n/sum(n))
  
  ready_for_overall$variable <- factor(ready_for_overall$variable,
                                       levels=rev(variables))
  
  # rename variables
  if(!is.null(variables_renames)) {
    ready_for_overall$variable <- recode(ready_for_overall$variable, !!!variables_renames)
  }
  
  # only get the classes of importance, it is assumed that it would have any problems since it was checked before
  ready_for_overall <- ready_for_overall[ready_for_overall$value %in% positive_values, ]
  
  # plot the overall population heatmap
  heatmap_overall <- ready_for_overall %>% 
    ggplot(aes(x=1, y=variable)) + 
    geom_tile(aes(fill=percent)) + 
    theme_minimal() + 
    scale_colour_stepsn(colors = palette, aesthetics = "fill",
                        limits=c(0, 1),
                        breaks=seq(0.1, 0.9, 0.1)) +
    ggtitle(' ') +
    theme(legend.position='right',
          #axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  # plot the combined heatmap (groups and population)
  heatmap_combined <- ggarrange(main_heatmap, heatmap_overall,
                                #labels=names(each_plot),
                                nrow=1, 
                                ncol=2,
                                widths=c(4, 1),
                                common.legend=F)
  
  if(!is.null(overall_heatmap_file)) {
    ggsave(output_overall_heatmap_file, plot=heatmap_combined)
  }
  
  ## heatmap R plot
  ro <- ready_for_overall[, c('variable', 'value', 'percent')] %>% as.data.frame
  ro$overall_percent <- ro$percent
  ro$percent <- NULL
  
  ready_both <- ready_for_heat %>%
    as.data.frame %>% 
    merge(ro,
          by=c('variable', 'value'),
          sort=F)
  
  ready_both$R <- ready_both$percent / ready_both$overall_percent
  
  # select only the data we want to plot
  ready_both <- ready_both[ready_both$value %in% positive_values, ]
  
  # steps
  step <- (max(ready_both$R) - min(ready_both$R))/10
  breaks_plot <- round(seq(min(ready_both$R) + step, max(ready_both$R) - step,  step), 2)
  
  heatmap_R <- ready_both %>% 
    ggplot(aes_string(x=predclass, y='variable')) + 
    geom_tile(aes(fill=R)) + 
    theme_minimal() + 
    scale_colour_stepsn(colors = palette, aesthetics = "fill", 
                        limits=c(min(ready_both$R), max(ready_both$R)),
                        breaks=breaks_plot) +
    ggtitle(' ') +
    theme(legend.position='right',
          #axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          #axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + theme(axis.text.y = element_text(hjust=0))
  
  if(!is.null(output_csv_values)) {
    write.csv(ready_both, output_csv_values, row.names=F)
  }
  if(!is.null(output_r_heatmap_file)) {
    ggsave(output_r_heatmap_file, plot=heatmap_R)
  }
  
  return(list(heatmap_r=heatmap_R,
              table_r=ready_both,
              heatmap_both=heatmap_combined,
              heatmap_main=main_heatmap,
              heatmap_overall=heatmap_overall))
}