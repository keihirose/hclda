hclda = function(lab, datx, r, N.press, type="exact", hierarchy=TRUE, trace=TRUE){
  if(class(lab)!="factor") stop('"lab" should be a "factor".')
  if(class(datx)[1]!="matrix") stop('"datx" should be a "matrix".')
  if(length(lab) != nrow(datx)) stop('The length of "lab" must be equal to the number of rows of "datx".')
  if(type!="exact" & type!="fast") stop('"type" should be "exact" or "fast".')
  p = ncol(datx)
  N = nrow(datx)
  
  lab.for = lab
  lab.latest = lab
  min_pair = NULL
  #min_err = LDA(lab, datx, r, N.press)$err
  G.prime = length(levels(lab))
  lab.latest_list = vector("list", G.prime-1)
  min_err_mtx = matrix(1,1,G.prime-1)   #変更点6/2
  
  Cluster = paste("c", 1:G.prime, sep="")

  #廣瀬追記
  if(type=="exact") SS.prime <- creasteSS_Rcpp(lab, datx)
  if(type=="fast") SS.prime <- creasteSS_fast_Rcpp(lab, datx, r)
  
  
      #１個抜きクロスバリデーション
        cvres = cv_Rcpp(datx, 
                lab,
                r, 
                N.press, 
                SS.prime,
                type)
        err = cvres$err
  
  lab.latest_list[[1]] = lab    
  min_err_mtx[1,1] = err/N    
  #---------------------------↑追記↑---------------------------#

  if(hierarchy){
    
    for(step in 1:(G.prime-2)){   #変更点6/2
      lab.step = lab.latest
      G.step = length(levels(lab.step))
      min_err = 1
      
      for(par1 in 1:(G.step-1)){
        for(par2 in (par1+1):G.step){
          lab.for = lab.step
          levels(lab.for)[c(par1,par2)] <- Cluster[step]

          #廣瀬追記
          if(type=="exact") SS.for <- creasteSS_Rcpp(lab.for, datx)
          if(type=="fast") SS.for <- creasteSS_fast_Rcpp(lab.for, datx, r)
          ## 一個抜きクロスバリデーション
          err = 0
          SS.for.cv <- vector(length(Cluster),mode="list")
          names(SS.for.cv) <- Cluster

         #クロスバリデーション！
         cvres = cv_clustering_Rcpp(datx, 
                  lab,
                  lab.for,
                  r, 
                  N.press, 
                  SS.for,
                  Cluster,
                  type)
          err = cvres$err


          #誤判別数を率に変換
          err_rate = err/N
          
          #誤判別率の更新とその時のペアを記録
          if(min_err > (err_rate)){
            min_err = err_rate
            min_pair = c(par1, par2)
            
            #ラベルを更新
            lab.latest = lab.step
            levels(lab.latest)[min_pair] = Cluster[step]
          }
        }
      }
      
      #クラスターが更新されなければ終了
      #if(G.step == length(levels(lab.latest))) break
      lab.latest_list[[step+1]] = lab.latest   #変更点6/2
      min_err_mtx[1,step+1] = min_err   #変更点6/2
      
      if(trace){
        # cat( "\014" )
        #cat(paste("d:", N.press,  "step:", step, "/"),file = "")
        cat(paste("step: ", step, "/", G.prime-2, ", ", sep=""),file = "")
        cat(c("cluster:", levels(lab.latest)),file = "")
        cat("\n")
      }
      
    }
  }
  
  #------------------↓三浦修正6/2↓------------------#
  min_num = which.min(min_err_mtx)
  best_error = min_err_mtx[min_num]
  best_cluster = lab.latest_list[[min_num]]
  
  best_cluster_lab = array(NA,G.prime)
  for(i in 1:G.prime){
    best_cluster_lab[i] = as.character((best_cluster[lab==levels(lab)[i]])[1])
  }
  best_cluster_f = factor(best_cluster_lab)
  names(best_cluster_f) = levels(lab)

    if(!hierarchy){
      min_err_mtx <- min_err_mtx[1]
      lab.latest_list <- lab.latest_list[[1]]
    }

  

  #廣瀬追加
  resultall <- list(err_rate=min_err_mtx, lab.latest=lab.latest_list, best_error=best_error, best_cluster=best_cluster, best_cluster_f=best_cluster_f,
              datx=datx, lab=lab, r=r, N.press=N.press) #r, N.pressを追加
  resultall$call <- match.call()
  class(resultall) <- "hclda"
  return(resultall)
  #------------------↑三浦修正6/2↑------------------#
}

