\name{pdf.figure.proof}
\alias{pdf.figure.proof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a nice pdf for a single plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pdf.figure.proof(findex = 1, file = paste("./pdf_sub_", findex, ".pdf", sep = ""), pw = 7, ph = 7, textscale = 0.6, xpd = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{findex}{
%%     ~~Describe \code{findex} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{pw}{
%%     ~~Describe \code{pw} here~~
}
  \item{ph}{
%%     ~~Describe \code{ph} here~~
}
  \item{textscale}{
%%     ~~Describe \code{textscale} here~~
}
  \item{xpd}{
%%     ~~Describe \code{xpd} here~~
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
function (findex = 1, file = paste("./pdf_sub_", findex, ".pdf", 
    sep = ""), pw = 7, ph = 7, textscale = 0.6, xpd = NA) 
{
    plotwidth <- pw/cm(1)
    plotheight <- ph/cm(1)
    margin <- 5/cm(1)
    pdfwidth <- plotwidth + 2 * margin
    pdfheight <- plotheight + 2 * margin
    posplot = c(margin/pdfwidth, (margin + plotwidth)/pdfwidth, 
        margin/pdfheight, (margin + plotheight)/pdfheight)
    pdf(file = file, height = pdfheight, width = pdfwidth, useDingbats = FALSE)
    par(cex = textscale, fig = posplot, mai = c(0, 0, 0, 0), 
        xpd = xpd)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
