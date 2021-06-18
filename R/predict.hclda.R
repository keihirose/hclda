


#loo_CLD_ez = function(test.dat, test.lab, train.dat, train.lab, r, N.press){
#  
#  n = nrow(train.dat)  #データのサンプル数
#  p = ncol(train.dat)  #ベクトルの次元数
#  G = length(levels(train.lab))  #群の数
#  
#  # n.g[g,1]でg群の要素数を表す 
#  n.g = matrix(summary(train.lab),G,1)
#  rownames(n.g) = levels(train.lab)
#  
#  # 平均ベクトル bar.x[,g]でg群の平均値 barxは全体の平均値
#  bar.x <- matrix(0,nrow = p,ncol = G)
#  for(i in 1:G){
#    bar.x[,i] = colMeans(train.dat[which(train.lab == levels(train.lab)[i]),])
#  }
#  colnames(bar.x) = levels(train.lab) 
#  barx <- colMeans(train.dat)
#  
#  # 分散共分散行列と群間・郡内変動行列
#  S.x <- vector("list",G)
#  S.x <- by(train.dat,train.lab,function(x) var(x))
#  #xjxjT <- by(train.dat,train.lab,function(x) as.matrix(x) %*% t(as.matrix(x)) ) 
#  #nminus1_S.x <- by(t(datx),label,function(x) var(x)) #廣瀬追加
#  
#  B <- matrix(0,p,p)
#  W <- matrix(0,p,p)
#  for(i in 1:G){
#    B = B + n.g[i,1]*(bar.x[,i] - barx) %x% t(bar.x[,i] - barx)
#    W = W + (n.g[i,1] - 1)*S.x[[i]]
#  }
#  
#  # 固有値・固有ベクトル
#  eigen_W <- eigen(W + r*diag(p), symmetric=TRUE)
#  Winv_half <- eigen_W$vector %*% diag(eigen_W$values^(-1/2)) %*% t(eigen_W$vector)
#  result_eigen <- eigen(Winv_half %*% B %*% Winv_half, symmetric=TRUE)
#  A <- Winv_half %*% result_eigen$vectors #射影ベクトル A[,d]でd番目の固有ベクトル
#  Val <- result_eigen$values
#  
#  # 射影先plot(t(Z))で圧縮してプロットできる
#  Z <- t(A[,1:N.press]) %*% t(test.dat - barx)
#  M <- t(A[,1:N.press]) %*% (bar.x - matrix(barx,p,G))
#  d <- colSums((matrix(Z,N.press,G) - M)^2)
#  hat.test.lab <- levels(train.lab)[which.min(d)]
#  
#  # 誤判別率 間違いであれば１増える
#  err = (hat.test.lab != test.lab)
#  
#  #出力
#  return(list(err=err, hat.lab=hat.test.lab))
#}





computeA = function(trainlab, trainx, r){
  G = length(levels(trainlab))
  p = ncol(trainx)
  n.g = matrix(summary(trainlab),G,1)
  rownames(n.g) = levels(trainlab)
  
  # 平均ベクトル bar.x[,g]でg群の平均値 barxは全体の平均値
  bar.x <- matrix(0,nrow = p,ncol = G)
  for(i in 1:G){
    bar.x[,i] = colMeans(trainx[which(trainlab == levels(trainlab)[i]),])
  }
  colnames(bar.x) = levels(trainlab) 
  barx <- colMeans(trainx)
  
  # 分散共分散行列と群間・郡内変動行列
  S.x <- vector("list",G)
  S.x <- by(trainx,trainlab,function(x) var(x))
  
  B <- matrix(0,p,p)
  W <- matrix(0,p,p)
  for(i in 1:G){
    B = B + n.g[i,1]*(bar.x[,i] - barx) %x% t(bar.x[,i] - barx)
    W = W + (n.g[i,1] - 1)*S.x[[i]]
  }
  
  # 固有値・固有ベクトル
  eigen_W <- eigen(W + r*diag(p), symmetric=TRUE)
  Winv_half <- eigen_W$vector %*% diag(eigen_W$values^(-1/2)) %*% t(eigen_W$vector)
  result_eigen <- eigen(Winv_half %*% B %*% Winv_half, symmetric=TRUE)
  A <- Winv_half %*% result_eigen$vectors #射影ベクトル A[,d]でd番目の固有ベクトル
  return(list(A = A, bar.x = bar.x))
}






predict.hclda = function(object, newx, ...){

  datx = object$datx
  lab = object$lab
  r = object$r
  N.press = object$N.press
  trainlab.cluster = object$best_cluster
  G.clu = length(levels(trainlab.cluster))  

  if(class(newx)[1]!="matrix") stop('"newx" should be a "matrix".')
  if(ncol(newx)!=ncol(datx)) stop('The dimension of new observations, "ncol(newx)", must be equal to that of the training data')
  # p = ncol(newx)
  N = nrow(newx)
  
  ### クラスターへの射影軸を計算 ###
  res = computeA(trainlab.cluster, datx, r)
  A = res$A
  bar.x = res$bar.x  
  
  ### クラスターごとの射影軸を計算 ###
  A.clu = list()
  bar.x.clu = list()
  for(clu in 1:G.clu){
    clustname = levels(trainlab.cluster)[clu]
    train.dat = datx[trainlab.cluster==clustname,]
    train.lab = factor(lab[trainlab.cluster==clustname])
  
    res.clu = computeA(train.lab, train.dat, r)
    A.clu[[clustname]] = res.clu$A
    bar.x.clu[[clustname]] = res.clu$bar.x
  }
  
  
  # テストデータとクラスター中心点の射影先
  Z <- newx %*% A[,1:N.press]
  M <- t(bar.x) %*% A[,1:N.press]
  
  # テストデータを一つずつ判別していく
  hat.testlab = array(NA,N)
  error = 0
  for(cv in 1:N){
    test.dat = as.matrix(t(newx[cv,]))
    Z.cv = Z[cv,]
    
    # 一回めの判別
    # result_clu = loo_CLD_ez(test.dat, test.lab, train.dat=datx, train.lab=trainlab.cluster, r, N.press)
    # hat.testlab[cv] = result_clu$hat.lab

    d <- rowSums((matrix(Z.cv, nrow(M), N.press, byrow=T) - M)^2)
    hat.testlab[cv] <- names(which.min(d))
    
    if(hat.testlab[cv] %in% trainlab.cluster){
      # 判別先がクラスターならばそのクラスター内部で再度判別
      # train.dat = datx[trainlab.cluster==result_clu$hat.lab,]
      # train.lab = factor(lab[trainlab.cluster==result_clu$hat.lab])
      # 
      # result_group = loo_CLD_ez(test.dat, test.lab, train.dat=train.dat, train.lab=train.lab, r, N.press=a)
      # error = error + result_group$err
      # hat.testlab[cv] = result_group$hat.lab
      
      Z.cv <- test.dat %*% A.clu[[hat.testlab[cv]]][,1:N.press]
      M.clu <- t(bar.x.clu[[hat.testlab[cv]]]) %*% A.clu[[hat.testlab[cv]]][,1:N.press]
      d <- rowSums((matrix(Z.cv, nrow(M.clu), N.press, byrow=T) - M.clu)^2)
      hat.testlab[cv] <- names(which.min(d))
    }
  }
  return(hat.lab=hat.testlab)
}
