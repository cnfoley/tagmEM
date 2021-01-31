# maximization step

m_step = function(data, clust_mn, clust_sd, clust_pi, alpha, K, junk_mixture = T, mu = 0){
  m = length(data);
  if(junk_mixture){
    K = K-1;
  }else{    
    mu = NULL;
  }
  den = colSums(alpha);
  if(K>1){
    tmp_mn = c(colSums(alpha[,1:K]*data)/den[1:K], mu);
    tmp_num = colSums(alpha[,1:K]*(matrix(rep(data,K) - rep(clust_mn[1:K],each=m), nrow = m))^2);
  }else if(K==1){
    tmp_mn = c(sum(alpha[,1:K]*data)/den[1:K], mu);
    tmp_num = sum(alpha[,1:K]*(matrix(rep(data,K) - rep(clust_mn[1:K],each=m), nrow = m))^2);
  }
  tmp_sig = sqrt(tmp_num/den[1:K]);
  tmp_pi = den/m;
  return(list("clust_mn" = tmp_mn, "clust_sd" = tmp_sig, "clust_pi" = tmp_pi))
}
