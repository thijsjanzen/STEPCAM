\name{plotSMC}
\alias{plotSMC}
\title{
  Plot the progression of the SMC algorithm.
}
\description{
  This function plots the progression over the several iterations of the ABC-SMC algorithm. It can be used to explore results obtained with the STEPCAM-ABC function.
}
\usage{
plotSMC(path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{Working Directory}
}

\author{Thijs Janzen}

\examples{
Artificial.Data <- generate.Artificial.Data(n_species = 40, n_traits = 3,
                        n_communities = 5, occurence_distribution = 0.2,
                        average_richness = 0.5, sd_richness = 0.2,
                        mechanism_random=FALSE)
O <- STEPCAM_ABC(Artificial.Data$abundances, Artificial.Data$traits,
                 numParticles = 10, n_traits = 3, plot_number = 1, stopRate = 0.8)
currentDir <- getwd();
plotSMC(paste(currentDir,"/",sep=""));

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SMC}
\keyword{ABC}
\keyword{community}
\keyword{assembly}% __ONLY ONE__ keyword per line
