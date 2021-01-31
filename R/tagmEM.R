#' Function to perform generalised t-augmented Gaussian mixture modelling using expectation-maximisation
#'
#' @param data numeric object of univariate observations
#' @param obs_names character vector of length the number of observations
#' @param max_iter numeric integer denoting the maximum number of iterations
#' before stopping the EM-algorithm's search for a maxima in the
#' log-likelihood.
#' @param tol numeric scalar denoting the maximum absolute difference between
#' two computations of the log-likelihood with which we accept that a maxima in
#' the log-likelihood has been computed.
#' @param cluster_sizes integer varying from 0 to number of data points N 
#' @param sigk_thresh lower bound of estimated cluster standard deviation. 
#' NOTE: avoids singular estimates of the 'sd'  
#' @param junk_mixture default is TRUE 
#' @param df default is 'df = 4'
#' @param junk_mean numeric scalar denoting the mean of the generalised
#' t-distribution. By default mean is set to zero.
#' @param junk_sd numeric scalar denoting the scale parameter in the generalised
#'  t-distribution
#' @param stop_bic_iter numeric integer I, for computational efficiency -
#' particularly when analysing large numbers of variants - we can stop the
#' EM-algorithm if the BIC is monotonic increasing over the previous I increases
#'  in the number of clusters K. By default evidence supporting at least 10
#'  clusters in the data is computed and so, for example, if the BIC from models
#'   which assume 6 clusters; 7 clusters; ... or; 10 clusters is monotonic
#'   increasing - in the number of clusters K -then the EM-algorithm is stopped
#'   and the model whose K minimises the BIC is returned.
#' @param min_clust_search numeric integer which denotes the minimum number of
#' clusters searched for in the data - default computes evidence supporting up
#' to K=10 clusters which might explain any clustered heterogeneity in the data.
#' @param results_list character list allowing users to choose whether to return
#'  a table with the variants assigned to: "all" of the clusters; a single
#'  "best" cluster or; both. By default we return both, i.e. results_list =
#'  list("all", "best").
#' @param  rand_sample random probability of being assignment to junk cluster   
#' @import dplyr reshape2 matrixStats data.table methods
#' @importFrom stats dnorm kmeans median sd
#' @return Returned are: estimates of the putative number of clusters in the
#' sample, allocation probabilities and summaries of the
#' association estimates for each observation;
#' @export 
tagmEM = function(data, obs_names = NULL, max_iter = 5e3, tol = 1e-5,
                   cluster_sizes = NULL, sigk_thresh = 1e-5,
                   junk_mixture = TRUE, df = 4, junk_mean = NULL, junk_sd = NULL,
                   stop_bic_iter = 5, min_clust_search = 8,
                   results_list = list("all", "best"),
                   rand_sample = seq(0.05,0.4, by = 0.05)){
  
  m = length(data);
  if(is.null(cluster_sizes)){cluster_sizes = 1:(min(9,m-1))}
  if(junk_mixture){
    cluster_sizes = c(0,cluster_sizes);
    if(is.null(junk_mean)){
      mu = median(data);
    }else{
      mu = junk_mean
    }
    if(is.null(junk_sd)){
      rng = range(data);
      sig = sd(data) + (rng[2]-rng[1]);
    }else{
      sig = junk_sd;
    }
  }
  num_clust = length(cluster_sizes);
  
  if (is.null(obs_names)) {
    obs_names <- paste0("ob_", 1:m)
  }
  
  
  bic_clust = vector();
  mn_clust = vector('list', num_clust);
  sd_clust = vector('list', num_clust);
  pi_clust = vector('list', num_clust);
  log_like = vector('list', num_clust);
  
  count0 = 1;
  
  for(K in cluster_sizes){
    if(K>0){
      kmn = kmeans(data,K);
      kmn_df = data.frame(x = data, cluster = kmn$cluster);
      kmn_sum_df = kmn_df %>% group_by(cluster) %>% 
        summarize(mn = mean(x), std = sd(x), size = n());
      
      kmn_sum_df = kmn_sum_df %>% mutate(alpha = size/sum(size));
      
      clust_mn = kmn_sum_df$mn;
      clust_sd = kmn_sum_df$std;
      clust_pi = kmn_sum_df$alpha;
      clust_sd[is.na(clust_sd)] = 1
      
      if(!junk_mixture){
        new_clust_mn = clust_mn;
        new_clust_sd = clust_sd;
        new_clust_pi = clust_pi;
        new_K = K; 
      }else{
        new_clust_mn = c(clust_mn, mu);
        new_clust_sd = clust_sd;
        junk_pr = sample(rand_sample, 1);
        new_clust_pi = c((1-junk_pr)*clust_pi, junk_pr);
        new_K = K+1;
      }
      count = 1;
      loglik_vec = NULL;
      loglik_diff = 1;
      while(count<max_iter & loglik_diff>tol){
        if(count == 1){
          estep = e_step(data, new_clust_mn, new_clust_sd, new_clust_pi, new_K, 
                          junk_mixture = junk_mixture, mu = mu, sig = sig, 
                          sigk_thresh = sigk_thresh);
          mstep = m_step(data, new_clust_mn, new_clust_sd, new_clust_pi, estep[["alpha"]], new_K,
                          junk_mixture = junk_mixture, mu = mu);
          loglik_vec = c(loglik_vec, estep[["loglik"]]);
        }else{
          loglik0 = estep[["loglik"]];
          estep = e_step(data, mstep[["clust_mn"]], mstep[["clust_sd"]], 
                          mstep[["clust_pi"]], new_K,
                          junk_mixture = junk_mixture, mu = mu, sig = sig,
                          sigk_thresh = sigk_thresh);
          mstep = m_step(data, mstep[["clust_mn"]], mstep[["clust_sd"]], 
                          mstep[["clust_pi"]], estep[["alpha"]], new_K,
                          junk_mixture = junk_mixture, mu = mu);
          if(is.na(estep[["loglik"]])){
            loglik_diff = tol/2;
          }else{
            loglik_diff = abs(loglik0 - estep[["loglik"]]); 
          }
          loglik_vec = c(loglik_vec, estep[["loglik"]]);
        }
        count = count +1;
      }
    }else{
      new_K = 1;
      mstep = estep = list();
      mstep[["clust_mn"]] = mu;
      mstep[["clust_sd"]] = NA;
      mstep[["clust_pi"]] = 1;
      estep[["loglik"]] = sum(gen_t_tagm(data, df = df, mu = mu, sig = sig, log=T));
    }
    mn_clust[[count0]] = mstep[["clust_mn"]];
    sd_clust[[count0]] = mstep[["clust_sd"]];
    pi_clust[[count0]] = mstep[["clust_pi"]];
    log_like[[count0]] = estep[["loglik"]];
    if(junk_mixture){
      bic_clust[count0] = -2*estep[["loglik"]] + (3*(new_K - 1))*log(m);
    }else{
      bic_clust[count0] = -2*estep[["loglik"]] + (3*new_K - 1)*log(m);
    }
    
    
    if((count0 > stop_bic_iter) & (K >= min_clust_search)){
      tmp_stop = which(!is.na(bic_clust));
      if(length(tmp_stop)>stop_bic_iter){
        tmp_cond = TRUE;
        rng_iter = tmp_stop[(length(tmp_stop)-stop_bic_iter):length(tmp_stop)];
        for(j in 1:stop_bic_iter){
          tmp_cond = tmp_cond & (bic_clust[rng_iter[j]] < bic_clust[rng_iter[j+1]]);
        }
        if(tmp_cond){
          break;
        }
      }
    }
    count0 = count0 +1;
  }
  
  results = clust_prob_tagm(bic_clust, data, mn_clust, sd_clust, pi_clust,
                            junk_mixture = junk_mixture, df = df, mu = mu, sig = sig, 
                            obs_names = obs_names,
                            sigk_thresh = sigk_thresh);
  results_all = results[order(results$cluster),];
  results_best = best_clust(results)
  res = list(results = list(all = results_all, best = results_best), 
             log_likelihood = log_like, bic = bic_clust)
  return(res)
}
