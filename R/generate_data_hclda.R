generate_data_hclda = function(N,bar.x_samp){
  
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