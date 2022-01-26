generate_data_hclda = function(N,bar.x_samp){

  if(missing(bar.x_samp)){
    p <- 20
    bar.x_samp <- cbind(c1 = c(rep(2,p/2), rep(0,p/2)),
                        c2 = c(rep(0,p/2), rep(2,p/2)),
                        c3 = c(rep(-1,p/2), rep(sqrt(3),p/2)),
                        c4 = c(rep(-1,p/2), rep(-sqrt(3),p/2)),
                        c5 = c(rep(0,p/2), rep(-2,p/2))
                        )
  }
  
  p = nrow(bar.x_samp)
  G = ncol(bar.x_samp)
  
  #ラベルの生成
  lab.prime = paste0("group",sort(ceiling(runif(n = N, min = 0, max = G))))
  lab = factor(lab.prime, levels = unique(lab.prime))
  #各群の個数
  n.g_samp = table(lab)
  
  #データの生成
  datx = matrix(0, nrow=N, ncol=p)
  n.g_ruikei = 0 #累計値
  for(g in 1:G){
    for(k in 1:p){
      datx[(n.g_ruikei+1):(n.g_samp[g]+n.g_ruikei),k] = rnorm(n=n.g_samp[g],bar.x_samp[k,g],1) 
    }
    n.g_ruikei = n.g_ruikei + n.g_samp[g]
  }
  #datx = as.data.frame(datx)
  
  
  return(list(lab=lab, datx=datx))
}