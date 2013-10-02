\name{PredictLogP}
\alias{PredictLogP}
\title{
  Predict the logP value of molecules in a given folder
}
\description{
  This function uses a pre-trained logP prediction model to estimate the logP value for a number of given molecules.
}
\usage{
  PredictLogP(structures.file, threads = -1)
}
\arguments{
  \item{structures.file}{A file containing the molecules. Can be SMILES or SDF.}
  \item{threads}{The number of threads to spawn when calculating the descriptors. Default is -1 and means use as many threads as are available.}
}
\details{
\code{PredictLogP} can be used to estimate the logP value of molecules. Descriptors are calculated using the Java based PaDEL-Descriptor package. Missing values are imputed using the \code{impute.knn} function in the \code{impute} package. These descriptors are used by a pre-trained support vector regression model to estimate logP values for new molecules. TBD: add details about the nature of the input format of the molecule files (see PD).
}
\value{
  A two column /code{data.frame} is returned with each row containing the name of the molecule as presented in the input file and the predicted logP value for that molecule.
}
\references{
TBD: insert reference to published paper
}
\author{
Daniel Murrell
}
\note{
TBD_note
}

\seealso{
TBD_seealso
%\code{\link{GenerateDescriptors}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%predictions <- PredictLogP(system.file("test_structures", "structures_10.sdf", package="camb"))
}
\keyword{ logP }
