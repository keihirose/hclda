\name{predict.hclda}
\alias{predict.hclda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the class using the hclda object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function gives us the prediction label for the new observations.   
}
\usage{
\method{predict}{hclda}(object, newx,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The hlcda object.
}
  \item{newx}{
Matrix of new observations.
}
  \item{\dots}{
Further arguments.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}



\value{
\item{hat.lab}{Prediction of labels for new observations}
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

\seealso{\code{hclda} objects}

\examples{
#generate data
set.seed(1)
dat <- generate_data_hclda(N=500) #generate dataset
lab <-  dat$lab
datx <-  as.matrix(dat$datx)

#fitting
fit <- hclda(lab, datx, r=1e-5, N.press=1, type="exact") #hclda with exact CV computation
fit.fast <- hclda(lab, datx, r=1e-5, N.press=1, type="fast") #hclda with approximate CV computation (fast computation)
fit
fit.fast

#prediction
newdata <- generate_data_hclda(N=100)
newx <- newdata$datx
pred <- predict(fit, newx)
pred
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
