\name{get_stats}
\alias{get_stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
for given clusters, outputs statistics of each cluster (number of cells, genes, average, top genes),
and outputs as a .csv matrix}
\usage{
get_stats(tenx, file_dest, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tenx}{
Input Seurat Object}
  \item{file_dest}{
destination of file}
  \item{d}{
name of file
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
function (tenx, file_dest, d)
{
    num_find <- 50
    aoe <- c("Group", "cell_number", "avg_gene", "avg_read",
        "avg_umi")
    for (i in 1:num_find) {
        aoe <- c(aoe, paste("top_", i, sep = ""))
    }
    df <- data.frame(aoe)
    for (groups in levels(tenx@ident)) {
        subgroup <- SubsetData(tenx, ident.use = groups)
        aod <- c(groups)
        aod <- c(aod, length(subgroup@data.info$nGene))
        aod <- c(aod, mean(subgroup@data.info$nGene))
        aod <- c(aod, sum(subgroup@raw.data)/length(subgroup@data.info$nGene))
        aod <- c(aod, mean(subgroup@data.info$nUMI))
        markers <- FindMarkers(tenx, groups)
        top_markers <- row.names(markers)[1:num_find]
        for (topm in top_markers) {
            aod <- c(aod, topm)
        }
        print(aod)
        df[groups] <- aod
        print(groups)
    }
    write.csv(df, paste(file_dest, d, "_stats.csv", sep = ""),
        row.names = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
