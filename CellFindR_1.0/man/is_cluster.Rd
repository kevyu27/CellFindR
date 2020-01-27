\name{is_cluster}
\alias{is_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cluster Determination 
}
\description{
asks if a particular clustering fulfill's CellFindR's metric:
the default is N-1 clusters fulfilling 10 genes with 2 fold expression}
\usage{
is_cluster(tenx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tenx}{
Seurat object}
}
\details{
this will base the clustering on the object@ident slicing of the metadata}