\name{summarise.glm}
\alias{summarise.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Writes a csv table summarising a list of models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summarise.glm(lstModels, outfunc = exp, writetab = TRUE, file = "modsum.csv", sigdigits = 3, transpose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lstModels}{
%%     ~~Describe \code{lstModels} here~~
}
  \item{outfunc}{
%%     ~~Describe \code{outfunc} here~~
}
  \item{writetab}{
%%     ~~Describe \code{writetab} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{sigdigits}{
%%     ~~Describe \code{sigdigits} here~~
}
  \item{transpose}{
%%     ~~Describe \code{transpose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lstModels, outfunc = exp, writetab = TRUE, file = "modsum.csv", 
    sigdigits = 3, transpose = FALSE) 
{
    nomods <- length(lstModels)
    allCoeffs <- c()
    for (i in 1:nomods) {
        mod <- lstModels[[i]]
        vars <- names(mod$coefficients)
        novars <- length(vars)
        for (j in 1:novars) {
            curname <- vars[j]
            var_present <- (curname \%in\% allCoeffs)
            if (!(var_present)) 
                allCoeffs <- c(allCoeffs, curname)
        }
    }
    noCoeffs <- length(allCoeffs)
    matPointEst <- matrix(NA, nrow = nomods, ncol = noCoeffs, 
        dimnames = list(1:nomods, allCoeffs))
    matLB <- matrix(NA, nrow = nomods, ncol = noCoeffs, dimnames = list(1:nomods, 
        allCoeffs))
    matUB <- matrix(NA, nrow = nomods, ncol = noCoeffs, dimnames = list(1:nomods, 
        allCoeffs))
    vecAIC <- vector(mode = "numeric", length = nomods)
    vecDEX <- vector(mode = "numeric", length = nomods)
    for (i in 1:nomods) {
        mod <- lstModels[[i]]
        cis <- confint.default(mod)
        vars <- names(mod$coefficients)
        novars <- length(vars)
        vecAIC[i] <- mod$aic
        vecDEX[i] <- (1 - mod$deviance/mod$null.deviance) * 100
        for (j in 1:novars) {
            curname <- vars[j]
            matPointEst[i, curname] <- mod$coefficients[curname]
            matLB[i, curname] <- cis[curname, 1]
            matUB[i, curname] <- cis[curname, 2]
        }
    }
    if (writetab) {
        if (transpose) {
            strTable <- ""
            strTable <- paste(strTable, "Parameter", sep = "")
            for (i in 1:noCoeffs) strTable <- paste(strTable, 
                ",", allCoeffs[i], ",", allCoeffs[i], sep = "")
            strTable <- paste(strTable, ",AIC,DEX\n", sep = "")
            strTable <- paste(strTable, "Model", sep = "")
            for (i in 1:noCoeffs) strTable <- paste(strTable, 
                ",PE,CI", sep = "")
            strTable <- paste(strTable, ",AIC,DEX\n", sep = "")
            for (i in 1:nomods) {
                strTable <- paste(strTable, i, sep = "")
                for (j in 1:noCoeffs) {
                  curname <- allCoeffs[j]
                  curPE <- signif(outfunc(matPointEst[i, curname]), 
                    digits = sigdigits)
                  curLB <- signif(outfunc(matLB[i, curname]), 
                    digits = sigdigits)
                  curUB <- signif(outfunc(matUB[i, curname]), 
                    digits = sigdigits)
                  if (is.na(curPE)) {
                    strTable <- paste(strTable, ",", "-", ",", 
                      "-", sep = "")
                  }
                  else {
                    strTable <- paste(strTable, ",", curPE, ",", 
                      "(", curLB, "--", curUB, ")", sep = "")
                  }
                }
                mod <- lstModels[[i]]
                curAIC <- round(mod$aic, digits = 1)
                curDEX <- round((1 - mod$deviance/mod$null.deviance) * 
                  100, digits = 1)
                strTable <- paste(strTable, ",", curAIC, ",", 
                  curDEX, "\n", sep = "")
            }
        }
        else {
            strTable <- ""
            strTable <- paste(strTable, ",Model 1", sep = "")
            if (nomods > 1) 
                for (i in 2:nomods) strTable <- paste(strTable, 
                  ",,Model ", i, sep = "")
            strTable <- paste(strTable, "\n", sep = "")
            if (nomods > 1) 
                for (i in 1:nomods) {
                  strTable <- paste(strTable, ",Estimate,(95\% CI)", 
                    sep = "")
                }
            strTable <- paste(strTable, "\n", sep = "")
            for (i in 1:noCoeffs) {
                curname <- allCoeffs[i]
                strTable <- paste(strTable, curname, sep = "")
                for (j in 1:nomods) {
                  curPE <- signif(outfunc(matPointEst[j, curname]), 
                    digits = sigdigits)
                  curLB <- signif(outfunc(matLB[j, curname]), 
                    digits = sigdigits)
                  curUB <- signif(outfunc(matUB[j, curname]), 
                    digits = sigdigits)
                  if (is.na(curPE)) {
                    strTable <- paste(strTable, ",", "-", ",", 
                      "-", sep = "")
                  }
                  else {
                    strTable <- paste(strTable, ",", curPE, ",", 
                      "(", curLB, "--", curUB, ")", sep = "")
                  }
                }
                strTable <- paste(strTable, "\n", sep = "")
            }
            strTable <- paste(strTable, "AIC", sep = "")
            for (i in 1:nomods) {
                mod <- lstModels[[i]]
                curAIC <- round(mod$aic, digits = 1)
                strTable <- paste(strTable, ",", curAIC, ",", 
                  sep = "")
            }
            strTable <- paste(strTable, "\n", sep = "")
            strTable <- paste(strTable, "DEX", sep = "")
            for (i in 1:nomods) {
                mod <- lstModels[[i]]
                curDEX <- round((1 - mod$deviance/mod$null.deviance) * 
                  100, digits = 1)
                strTable <- paste(strTable, ",", curDEX, ",", 
                  sep = "")
            }
            strTable <- paste(strTable, "\n", sep = "")
        }
        cat(strTable, file = file)
    }
    data.frame(pe = matPointEst, lb = matLB, ub = matUB, aic = vecAIC, 
        dex = vecDEX)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
