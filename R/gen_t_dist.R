# generalised t-dist

gen_t_tagm = function(x, df = 4, mu, sig, log = FALSE){
  tmp = (gamma((df+1)/2)/gamma(df/2)/sqrt(pi*df)/sig)*(1+(x-mu)^2/sig^2/df)^(-(df+1)/2);
  if(log){tmp = log(tmp);}
  return(tmp);
}
