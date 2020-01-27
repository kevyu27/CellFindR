\name{plottingj}
\alias{plottingj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting genes
}
\description{
plots top 100 genes in violin and cluster plot }
\usage{
plottingj(sub_obj, file_dest, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sub_obj}{
Seurat Object}
  \item{file_dest}{
find destination}
  \item{d}{
name of file output}
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
function (sub_obj, file_dest, d)
{
    file_create <- paste(paste(file_dest, d, sep = ""))
    markers <- FindAllMarkers(sub_obj, only.pos = TRUE, min.pct = 0.25,
        thresh.use = 0.1)
    markers <- markers[order(markers$avg_diff, decreasing = TRUE),
        ]
    markers <- na.omit(markers[1:100, ])
    a <- (markers[, 6])
    dir.create(file_create)
    pdf(paste(file_create, "/TSNE.pdf", sep = ""))
    TSNEPlot(sub_obj)
    dev.off()
    dir.create(paste(file_create, "Cluster", sep = "/"))
    dir.create(paste(file_create, "Violin", sep = "/"))
    for (i in a) {
        pdf(paste(paste(file_create, "Cluster", "", sep = "/"),
            i, ".pdf", sep = ""))
        FeaturePlot(sub_obj, c(i), cols.use = c("grey", "blue"),
            pt.size = 2)
        dev.off()
        pdf(paste(paste(file_create, "Violin", "", sep = "/"),
            i, ".pdf", sep = ""))
        VlnPlot(sub_obj, c(i))
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
