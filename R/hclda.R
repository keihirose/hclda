parhclda <- function(npar, 
                     parALL, 
                     lab.step,
                     Cluster,
                     step,
                     type,
                     datx,
                     r,
                     lab,
                     N.press,
                     N){
  par1 <- parALL[npar,1]
  par2 <- parALL[npar,2]
  lab.for = lab.step
  levels(lab.for)[c(par1,par2)] <- Cluster[step]

  #廣瀬追記
  if(type=="exact") SS.for <- creasteSS_Rcpp(lab.for, datx)
  if(type=="fast") SS.for <- creasteSS_fast_Rcpp(lab.for, datx, r)
  ## 一個抜きクロスバリデーション
  #err = 0
  #SS.for.cv <- vector(length(Cluster),mode="list")
  #names(SS.for.cv) <- Cluster

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
  err_rate
}

hclda = function(lab, datx, r, N.press, type="exact", hierarchy="cv", parallel=FALSE, max.cores=10, trace=TRUE){
  if(class(lab)!="factor") stop('"lab" should be a "factor".')
  if(class(datx)[1]!="matrix") stop('"datx" should be a "matrix".')
  if(length(lab) != nrow(datx)) stop('The length of "lab" must be equal to the number of rows of "datx".')
  if(type!="exact" & type!="fast") stop('"type" should be "exact" or "fast".')
  if(hierarchy!="cv" & hierarchy!="ward" & hierarchy!="none") stop('"hierarchy" should be "cv", "none" or "ward".')

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

  if(hierarchy=="cv"){
    
    for(step in 1:(G.prime-2)){   #変更点6/2
      lab.step = lab.latest
      G.step = length(levels(lab.step))
      min_err = 1
      
      if(!parallel){
        for(par1 in 1:(G.step-1)){
          for(par2 in (par1+1):G.step){
            lab.for = lab.step
            levels(lab.for)[c(par1,par2)] <- Cluster[step]

            #廣瀬追記
            if(type=="exact") SS.for <- creasteSS_Rcpp(lab.for, datx)
            if(type=="fast") SS.for <- creasteSS_fast_Rcpp(lab.for, datx, r)
            ## 一個抜きクロスバリデーション
            #err = 0
            #SS.for.cv <- vector(length(Cluster),mode="list")
            #names(SS.for.cv) <- Cluster

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
        
        ##クラスターが更新されなければ終了
        ##if(G.step == length(levels(lab.latest))) break
        #lab.latest_list[[step+1]] = lab.latest   #変更点6/2
        #min_err_mtx[1,step+1] = min_err   #変更点6/2
        #
        #if(trace){
        #  # cat( "\014" )
        #  #cat(paste("d:", N.press,  "step:", step, "/"),file = "")
        #  cat(paste("step: ", step, "/", G.prime-2, ", ", sep=""),file = "")
        #  cat(c("cluster:", levels(lab.latest)),file = "")
        #  cat("\n")
        #}
        
      }


      if(parallel){
        parALL <- matrix(NA,G.step*(G.step-1)/2,2)
        count <- 0
        for(par1 in 1:(G.step-1)){
          for(par2 in (par1+1):G.step){
            count <- count + 1
            parALL[count,] <- c(par1,par2)
          }
        }
        nparALL <- count
        ncores <- min(detectCores(), max.cores, nparALL)
        scl = makeCluster(getOption("mc.cores", ncores), type = "FORK")
        #export_name <- c("parALL",
        #                 #"creasteSS_Rcpp", 
        #                 #"creasteSS_fast_Rcpp", 
        #                 #"cv_clustering_Rcpp",
        #                 "lab.step",
        #                 "Cluster",
        #                 "step",
        #                 "type",
        #                 "datx",
        #                 "r",
        #                 "lab",
        #                 "N.press",
        #                 "N")
        #clusterExport(scl, export_name, envir=environment())
        #clusterEvalQ(scl, library(vrec))
        #registerDoParallel(scl)
        result_tmp = clusterApplyLB(cl = scl, x=seq(1, nparALL), fun=parhclda, parALL=parALL,
                                                                               lab.step=lab.step,
                                                                               Cluster=Cluster,
                                                                               step=step,
                                                                               type=type,
                                                                               datx=datx,
                                                                               r=r,
                                                                               lab=lab,
                                                                               N.press=N.press,
                                                                               N=N)
        stopCluster(scl)
        min_err = min(unlist(result_tmp))
        min_pair = parALL[which.min(result_tmp), ]
        #ラベルを更新
        lab.latest = lab.step
        levels(lab.latest)[min_pair] = Cluster[step]
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
  
  if(hierarchy=="ward"){
    ## Ward's method ####
    means = apply(datx,2,function(x) tapply(x,lab,mean))
    dm = dist(means)
    hc = hclust(d=dm,method="ward.D2")

    for(step in 1:(G.prime-2)){
      Cluster.ward = cutree(tree=hc,k=G.prime-step)
      nclusters <- max(Cluster.ward)
      index_duplicated <- rep(FALSE, nclusters)
      for(ii in 1:nclusters) index_duplicated[ii] <- sum(Cluster.ward == ii) > 1
      n_duplicated <- sum(index_duplicated)
      index_tmp <- vector(n_duplicated, mode="list")
      which_index_duplicated <- which(index_duplicated)
      for(ii in 1:n_duplicated) index_tmp[[ii]] <- Cluster.ward == which_index_duplicated[ii]

      lab.for = lab
      #for(ii in 1:n_duplicated) levels(lab.for)[index_tmp[[ii]]] = Cluster[ii]

      for(ii in 1:n_duplicated){
        index_tmp2 <- levels(lab.for) %in% names(index_tmp[[ii]])[index_tmp[[ii]]]
        levels(lab.for)[index_tmp2] <- Cluster[ii]
      }


      #廣瀬追記
      if(type=="exact") SS.for <- creasteSS_Rcpp(lab.for, datx)
      if(type=="fast") SS.for <- creasteSS_fast_Rcpp(lab.for, datx, r)
      ## 一個抜きクロスバリデーション
      #err = 0
      #SS.for.cv <- vector(length(Cluster),mode="list")
      #names(SS.for.cv) <- Cluster

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

      lab.latest_list[[step+1]] = lab.for   #変更点1/12
      min_err_mtx[1,step+1] = err_rate   #変更点1/12

      if(trace){
        # cat( "\014" )
        #cat(paste("d:", N.press,  "step:", step, "/"),file = "")
        cat(paste("step: ", step, "/", G.prime-2, ", ", sep=""),file = "")
        cat(c("cluster:", levels(lab.for)),file = "")
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

    if(hierarchy=="none"){
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

