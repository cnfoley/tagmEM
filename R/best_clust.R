# function to identify the 'best' clustering of the data using BIC criterion

best_clust = function(dta){
  res = NULL;
  unique_obs = unique(dta$observation);
  for(i in unique_obs){
    tmp = dta$observation==i;
    mx = which.max(dta$probability[tmp]);
    res = rbind(res, dta[tmp,][mx,])
  }
  return(res);  
}

# prevents R check throwing errors about unassigned global variables
utils::globalVariables(c("cluster", "size", "x"))