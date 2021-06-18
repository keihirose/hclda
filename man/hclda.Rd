\name{hclda}
\alias{hclda}
\alias{print.hclda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical clustered linear discriminant analysis.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function gives us hierarchical clustered linear discriminant analysis.  
It can accurately classify data with a large number of groups, which is difficult to be classified with an ordinary linear discriminant analysis.
}
\usage{
hclda(lab, datx, r, N.press, type = "exact", hierarchy = TRUE, trace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lab}{
  Group label for datx (e.g., the name of the group.) The class of the lab object should be \code{factor}.
  }
  \item{datx}{
  Input matrix with each row being an observation vector. The number of the rows of "datx" must be equal to the length of "lab".
  }
  \item{r}{
  Regularization parameters for the within-class covariance matrix.
  }
  \item{N.press}{
  The number of dimension of the projected space when applying linear discriminant analysis.
  }
  \item{type}{
  Either \code{"exact"} or \code{"fast"}. If \code{"exact"}, the exact cross-validation value is computed.  The exact calculation becomes slow when the number of dimension of the observation is large.  In such cases, \code{"fast"} is much faster but the approximation is made.  The approximation error is small when the number of observations is sufficiently large.
  }
  \item{hierarchy}{
  If \code{TRUE}, the hierarchical clustering is conducted.  Otherwise, the ordinary LDA with cross-validation is performed.
  }
  \item{trace}{
  If \code{TRUE}, tracing information is printed.  Default is \code{TRUE}.
  }
%  \item{\dots}{Other graphical parameters to plot}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
  \item{err_rate}{
  Misclassification rate per step calculated by cross-validation. 
  }
  \item{lab.latest}{
  Cluster labels for each step.
  }
  \item{best_error}{
  Misclassification rate for the best step.
  }
  \item{best_cluster}{
  Cluster labels for the best step.
  }
  \item{best_cluster_f}{
  Allocation from groups to clusters for the best step.
  }
  \item{datx}{
  Input matrix.
  }
  \item{lab}{
  Group label for datx.
  }
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Kanta Miura <\email{kantamiura903@gmail.com}> and Kei Hirose <\email{mail@keihirose.com}>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{predict.hclda} object}

\examples{
#generate data
set.seed(1)
p = 20
bar.x_samp =cbind(
  c1 = c(rep(2,p/2), rep(0,p/2)),
  c2 = c(rep(0,p/2), rep(2,p/2)),
  c3 = c(rep(-1,p/2), rep(sqrt(3),p/2)),
  c4 = c(rep(-1,p/2), rep(-sqrt(3),p/2)),
  c5 = c(rep(0,p/2), rep(-2,p/2))
)
dat <- generate_data_hclda(N=500, bar.x_samp) #generate dataset
lab <-  dat$lab
datx <-  as.matrix(dat$datx)

#fitting
fit <- hclda(lab, datx, r=1e-5, N.press=1, type="exact") #hclda with exact CV computation
fit.fast <- hclda(lab, datx, r=1e-5, N.press=1, type="fast") #hclda with approximate CV computation (fast computation)
fit
fit.fast

#prediction
newdata <- generate_data_hclda(N=100, bar.x_samp)
newx <- newdata$datx
pred <- predict(fit, newx)
pred
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
