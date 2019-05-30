\name{getmatrix}
\alias{getmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generates differential gene matrix}
\description{
generates matrix of differential expression for given Seurat object with defined clusters: columns for each cluster: mean, std, avg diff, p_value (adjusted)}
\usage{
getmatrix(tenx, file_dest, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tenx}{
Seurat object
}
  \item{file_dest}{
location of destination file output
}
  \item{d}{
name of file }
}
\details{
The output files: will have in the columns for each cluster
mean, standard deviation of the mean, average differentation and p-value, all in the log space. }