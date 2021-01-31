# function to return cluster summaries for each number of clusters in the data K

clust_prob_tagm = function(bic_clust, data, mn_clust, sd_clust, pi_clust,
                           junk_mixture = TRUE, df = 4, mu = 0, sig = 1, 
                           obs_names, sigk_thresh){
  # Identify NAs from cluster parameter data
  keep_iter = NULL;
  for(j in 1:length(bic_clust)){
    tmp_keep = sum(is.na(mn_clust[[j]])) + sum(is.na(sd_clust[[j]])) +
      sum(is.na(pi_clust[[j]]));
    if(tmp_keep==0){
      keep_iter = c(keep_iter, TRUE);
    }else{
      keep_iter = c(keep_iter, FALSE);
    }
  }
  mx = which(bic_clust==min(bic_clust[which(!is.na(bic_clust)  & keep_iter)]));
  theta_mn = mn_clust[[mx]];
  clust_probs = pi_clust[[mx]];
  clust_sigs = sd_clust[[mx]];
  n_clust = length(theta_mn);
  m = length(data);
  e_step = e_step(data, theta_mn, clust_sigs, clust_probs, n_clust, 
                  junk_mixture = junk_mixture, mu = mu, sig = sig,
                  sigk_thresh = sigk_thresh);
  res1  = round(as.numeric(t(e_step[["alpha"]])),3);  
  res2 = round(rep(theta_mn, m),3);
  if(junk_mixture){
    clust_sigs = c(clust_sigs,NA)
  }
  res = data.table(rep( obs_names, each = n_clust), rep(1:n_clust, m), res1, res2, 
                   rep( data, each = n_clust), rep( clust_sigs, m))
  names(res) = c("observation", "cluster", "probability", "cluster_mean", "data", "clust_sd");
  factor_clust = res$cluster;
  if(junk_mixture){
    tmp = res$cluster == n_clust;
    factor_clust[tmp] = "Junk";
  }
  res$cluster_class = as.character(factor_clust);
  return(res)
}
