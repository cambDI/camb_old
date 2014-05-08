\name{StandardiseMolecules}
\alias{StandardiseMolecules}
\title{
Read and standardise molecules from an input SDF or SMILES structures file.
}
\description{
Reads and applies basic filtering operations to molecules from an input SDF or SMILES structures file. These include removal of non fully organic molecules, limiting the number of specific halogens, limiting the molecular mass of a molecule, striping out salts, standardising aromatic description and representing all tautomers in a single tautomeric form.  
}
\usage{
StandardiseMolecules(structures.file, standardised.file, removed.file = "", output = "standardisation_info.csv", remove.inorganic = FALSE, fluorine.limit = -1, chlorine.limit = -1, bromine.limit = -1, iodine.limit = -1, min.mass.limit = -1, max.mass.limit = -1, number.processed = -1)
}
\arguments{
  \item{structures.file}{
the name of the file in which the structures reside. Can be SDF or SMILES.
}
  \item{standardised.file}{
the name of the file which holds the standardised structures. Will be saved in SDF format.
}
  \item{removed.file}{
the name of the file where the structures of the molecules removed by the filtering parameters are stored.
}
  \item{output}{
the name of the file that stores additional property information that was included in the original structures file (SDF molecule attributes for example).
}
  \item{remove.inorganic}{
whether or not molecule which contain atoms not in {H, C, N, O, P, S, F, Cl, Br, I} should be removed.
}
  \item{fluorine.limit}{
if a molecule contains more than \code{fluorine.limit} atoms, that molecule is removed
}
  \item{chlorine.limit}{
if a molecule contains more than \code{chlorine.limit} atoms, that molecule is removed
}
  \item{bromine.limit}{
if a molecule contains more than \code{bromine.limit} atoms, that molecule is removed
}
  \item{iodine.limit}{
if a molecule contains more than \code{iodine.limit} atoms, that molecule is removed
}
  \item{min.mass.limit}{
if a molecule has a molecular mass less than \code{min.mass.limit}, that molecule is removed
}
  \item{max.mass.limit}{
if a molecule has a molecular mass greater than \code{max.mass.limit}, that molecule is removed
}
  \item{number.processed}{
the number of molecule to be processed by the standardisation routine. This is used to test the function on a few molecule in large structure files.
}
}
\details{
TBD
}
\value{
TBD something about the value returned
}
\references{
TBD put some references
}
\author{
Daniel Murrell, dsmurrell@gmail.com and Isidro Cortes <isidrolauscher@gmail.com>
}
\note{
TBD do we need notes here
}

\seealso{
TBD link to other code here. %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
TBD put example
}
\keyword{ ~StandardiseMolecules }