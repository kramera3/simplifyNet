% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GlobalSparse.R
\name{gns}
\alias{gns}
\title{Global Threshold Sparsification}
\usage{
gns(E_List, remove.prop, cutoff)
}
\arguments{
\item{E_List}{Edge list of the given network in the format of | node 1 | node 2 | weight |.}

\item{remove.prop}{The proportion of highest weighted edges to retain. A value between 0 and 1.}

\item{cutoff}{Threshold value for edge weight thresholding.}
}
\value{
Edge list of sparsified network
}
\description{
Remove all edges under certain edge weight threshold.
}
\examples{
g = igraph::erdos.renyi.game(100, 0.1)
igraph::E(g)$weight <- runif(length(igraph::E(g)))
E_List = cbind(igraph::as_edgelist(g), igraph::E(g)$weight)
H = gns(E_List, remove.prop = 0.5)
}
\author{
Alexander Mercier
}
