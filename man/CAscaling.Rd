\name{CAscaling}
\alias{CAscaling}

\title{
Performs CA scaling.
}

\description{CA scaling checks which item pairs are locally dependent. (Straat, Van der Ark, and Sijtsma, 
submitted). } 

\usage{ CAscaling(Data, Range=3) }

\arguments{
  \item{Data}{Matrix of integers, missing values are not allowed}
  \item{Range}{Range of the boxplot}
}

\value{
  \item{Output}{A matrix consisting of ncol(Data) columns and 4 rows. 
  Each column refers to an item.
  The first three rows summarize the results for inter-item covariances conditional on single item scores (row 1),
  inter-item covariances conditional on the rest score to detect negative local dependence (row 2),
  and inter-item covariances conditional on the rest score to detect positive local dependence (row 3).
  The fourth row shows which items are scalable (0) and unscalable (>0).}
  
  \item{Scores}{Score distributions corresponding to each of the first three rows in Output.}
  
  \item{Fences}{Upper fences of the boxplot, which are used to determine which item pairs are locally dependent.}
}


\author{
 J. H. Straat \email{j.h.straat@uvt.nl}
}

\references{
 Straat J. H., Van der Ark, L. A., and Sijtsma K. (2012). 
 Item scaling using conditional association.
 Manuscript submitted for publication.
}

\examples{
# Data
X <- matrix(sample(0:1,2000,TRUE),200,10)
dimnames(X)[[2]] <- paste("V",1:10,sep="")

# CA scaling
CAscaling(X)[[1]]
}
\keyword{CAscaling}
