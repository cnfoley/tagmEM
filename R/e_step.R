# expectation step

e_step = function(data, clust_mn, clust_sd, clust_pi, K, 
                  junk_mixture = T, mu = 0, sig = 1, sigk_thresh){
  m = length(data);
  if(junk_mixture){
    K = K-1;
  }
  tmp_theta = rep(data,K);
  tmp_mn = rep(clust_mn[1:K], each = m);
  tmp_sd = rep(clust_sd, each = m);
  tmp_sd[tmp_sd<sigk_thresh] = sigk_thresh;
  tmp_pi = rep(clust_pi[1:K], each = m);
  rjk_num = dnorm(tmp_theta, tmp_mn, tmp_sd, log=T) + log(tmp_pi);
  if(junk_mixture){
    rjk_num = c(rjk_num, gen_t_tagm(data, mu = mu ,sig = sig,log=T) + log(clust_pi[K+1]));
  }
  rjk_num_mat = matrix(rjk_num, nrow = m);
  rjk_row_mx = rowMaxs(rjk_num_mat);
  rjk_num_mat_std = rjk_num_mat - rjk_row_mx;
  rjk = exp(rjk_num_mat_std - log(rowSums(exp(rjk_num_mat_std))));
  log_lik = sum(rjk_row_mx + log(rowSums(exp(rjk_num_mat_std))));
  
  return(list("alpha" = rjk, "loglik" = log_lik))
}
