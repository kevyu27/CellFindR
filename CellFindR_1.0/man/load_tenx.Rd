\name{load_tenx}
\alias{load_tenx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
loads 10x matrix files and sets up as a Seurat object
}
\description{
takes in 10x matrix files and creates Seurat object in the same directory as variable name tenx
}
\usage{
load_tenx(file_loc, res = 1, proj_name = "tenx_data", cutoff = 7000, mito = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_loc}{
source directory folder holding matrix, barcodes and gene names files.
}
  \item{res}{
set initial resolution of clustering
}
  \item{proj_name}{
set name of project
}
  \item{cutoff}{
set maximum number of genes a cell can have, (doublet cut off)
}

  \item{mito}{
cut off 5 percent mito (True or false), default is false
}
}
