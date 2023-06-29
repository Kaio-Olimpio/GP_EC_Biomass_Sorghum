## Function recod
##
##' This function performs the recodification of a specific factor 
##' in the data set. It adds a new column with the coded levels to
##' the data set
##'
##' @title Easy factor recodification
##' @param data dataframe A dataframe containing the factor to be recoded
##' @param name.fact string Name of the column containing the factor to be recoded
##' @param cod string String that will appear before the numbers in the new code. Ex: "G001"
##' @return The function returns the datagrame with a new column containing the recoded levels
##' 
##' @author Saulo F. S. Chaves (saulo.chaves at ufv.br)
##' 


recod = function(data, name.fact, cod){
  stopifnot("cod must be a string" = is.character(cod))
  stopifnot("name.fact must be a string" = is.character(name.fact))
  
  num = dim(as.matrix(unique(data[,name.fact])))[1]
  data$ord = 1:nrow(data)
  
  a = data.frame(
    "JJJJ" = unique(data[,name.fact]),
    'XYWFX' = paste0(cod,sprintf(paste0('%0', nchar(num),'d'), seq(1:num)))
  ) 
  
  colnames(a)[which(colnames(a) == 'JJJJ')] = name.fact
  
  data = merge(data, a, by = name.fact, sort=F)
  colnames(data)[which(colnames(data) == "XYWFX")] = paste("cod",name.fact,sep = '.')
  
  data = data[order(data$ord),]
  data = data[,-which(colnames(data)=='ord')]
  rownames(data) = NULL
  
  return(data)
  
}
