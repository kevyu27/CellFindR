\name{just_plot}
\alias{just_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plots TSNE density plot of given gene
}
\description{
for input gene, looks at expression plotting using TSNE plot
}
\usage{
just_plot(tenx, gene, size = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tenx}{
input Seurat object
}
  \item{gene}{
name of gene
}
  \item{size}{
size of plotting dot
}
}
\details{
short-cut for the FeaturePlot from seurat}