\name{STEPCAM_ABC}
\alias{STEPCAM_ABC}
\title{
  ABC-SMC inference of the STEPCAM model
}
\description{
 ABC-SMC inference of the relative contribution of stochasticy (or random removal), filtering and limiting similarity (e.g. of shared resources (competition) or predators) through the use of stepwise commmunity assembly and fitting Functional Diversity patterns of simulated communities with those of observed communities. It takes quite a lot of computing time to run this function, so it is recommended to first run it for a very simple (species poor, artificial) dataset (see example).
}
\usage{
STEPCAM_ABC(data_abundances,data_species,numParticles,n_traits,plot_number,stopRate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_abundances}{
  \code{data_abundances} is a matrix containing the abundances of all species over all plots. \cr
  \code{rows}{ every row corresponds to a different plot or community} \cr
  \code{columns}{ every column corresponds to a different species.} \cr
  For an example of a suitable \code{data_abundances} matrix, run the \code{generate.Artificial.Data} function
}
  \item{data_species}{
  \code{data_species} is a matrix containing the species' specific traits. \cr
  \code{rows}{ every row corresponds to a different species} \cr
  \code{columns}{ the first column contains the name of the species, every consecutive column corresponds to a different trait} \cr
  It is important to realize that ideally, the number of traits is lower than the number of species in the most species poor plot of the 'data_abundances' matrix, as FRic calculation requires at least as many species as their are traits, and will autmatically omit some traits if this assumption is not met. \cr
  For an example of a suitable data_species matrix, run the generate.Artificial.Data function
}
  \item{numParticles}{The number of particles to be used in the ABC-SMC analysis.
}
  \item{n_traits}{The number of traits used in the analysis.
}
  \item{plot_number}{Row number of the plot for which we want to infer the contributions of stochasticity, filtering and competition.
}
  \item{stopRate}{The SMC procedure stops as soon as the acceptance rate of newly proposed parameter combinations drops below this level. Beware that low values
                  for this parameter might lead into long computation times.}
}

\value{
    \item{DA}{Posterior distribution of the relative contribution of dispersal assembly (random removal steps) to community assembly}
    \item{HF}{Posterior distribution of the relative contribution of habitat filtering to community assembly}
    \item{LS}{Posterior distribution of the relative contribution of limiting similarity to community assembly}
}
\references{ Fons van der Plas, Thijs Janzen, Alejandro Ordonez, Wimke Fokkema, Josephine Reinders, Rampal S. Etienne, and Han Olff 2015. A new modeling approach estimates the relative importance of different community assembly processes. Ecology 96:1502-1515.
}
\author{
Thijs Janzen & Fons van der Plas + Torsten Hauffe

}


\examples{
Artificial.Data <- generate.Artificial.Data(n_species = 40, n_traits = 3,
                        n_communities = 5, occurence_distribution = 0.2,
                        average_richness = 0.5, sd_richness = 0.2,
                        mechanism_random=FALSE)
O <- STEPCAM_ABC(Artificial.Data$abundances, Artificial.Data$traits,
                 numParticles = 10, n_traits = 3, plot_number = 1, stopRate = 0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ community }
\keyword{ assembly }
\keyword{ traits }
\keyword{ simulation }
\keyword{ processes }
\keyword{ functional }
\keyword{ diversity }% __ONLY ONE__ keyword per line
