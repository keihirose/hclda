print.hclda=function(x, digits = max(3, getOption("digits") - 3),num.result = 20,...){
		#check digits
		if(mode(digits)!="numeric")	stop('"digits" must be numeric.')
		if(length(digits) > 1)	stop('"digits" must be a scalar (1-dimensional vector).')
		if(as.integer(digits)!=digits)	stop('"digits" must be integer.')
		if(digits <= 0)	stop('"digits" must be positive integer.')

		#check num.result
		if(mode(num.result)!="numeric")	stop('"num.result" must be numeric.')
		if(length(num.result) > 1)	stop('"num.result" must be a scalar (1-dimensional vector).')
		if(as.integer(num.result)!=num.result)	stop('"num.result" must be integer.')
		if(num.result <= 0)	stop('"num.result" must be positive integer.')

		
		 err_rate <- x$err_rate
		 best_cluster <- x$best_cluster
		 best_cluster_f <- x$best_cluster_f



#RESULT
 	cat("\nCall:", paste(deparse(x$call)), "\n")
    cat("\nError rate:\n"); print(err_rate,digits=digits);
    cat("\nBest cluster:\n"); print(best_cluster_f,digits=digits);
   # cat("\nloss:", paste(x$loss));
    cat("\n")
 	invisible(x)      
    
}
