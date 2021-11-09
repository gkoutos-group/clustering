if(!("ggbiplot" %in% (.packages()))){
    library(ggbiplot)
}

pca_for_variables <- function(dataset, variables, target) {
    noNA <- rowSums(is.na(dataset[, variables])) == 0 #filter rows with missing values
    c_target <- target[noNA]
    dataset <- as.data.frame(lapply(dataset[noNA, c(variables)], as.numeric)) #get the numeric types of the columns
    
    pca <- prcomp(dataset[, variables], 
              center=T, 
              scale.=T)
    
    p <- ggbiplot::ggbiplot(pca, 
                             obs.scale=1, 
                             var.scale=1, 
                             groups=c_target, 
                             ellipse=T, 
                             circle=T) + 
    scale_color_discrete(name = '') + 
    theme(legend.direction = 'horizontal', legend.position='top')
    
    return(list(pca=pca, plot=p))
}