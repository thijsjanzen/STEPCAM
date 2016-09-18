\name{generate.Artificial.Data}
\alias{generate.Artificial.Data}
\title{
  Generate Artifical data that can be used by the STEPCAM model
}
\description{
  This function generates artifical data with a given number of species, traits and communities. The generated data can directly be used in
  STEPCAM-ABC inference.
}
\usage{
generate.Artificial.Data(n_species, n_traits, n_communities,
occurence_distribution, average_richness, sd_richness, mechanism_random)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_species}{The number of species in the species pool (so across all communities) of the desired dataset. }
  \item{n_traits}{The number of traits in the desired dataset.}
  \item{n_communities}{The number of communities in the desired dataset.}
  \item{occurence_distribution}{species occurrence distribution (analogous to species abundance distribution) in species pool. A lognormal distribution is assumed. The higher the value, the more uneven the occurrence distribution in the species pool}
  \item{average_richness}{average richness in local community as a proportion of richness in species pool}
  \item{sd_richness}{standard deviation in richness across communities}
  \item{mechanism_random}{if TRUE, species are randomly drawn from from species pool. If FALSE, it is assumed that a certain trait combination is optimal in each community and species with traits similar to this optimum more likely end up in the local community. Therefore, random.mechanism = FALSE assumes filtering processes to be most dominant in community assembly}
}
\value{
  \item{traits}{A matrix containing the traits for all species, rows are the different species, column correspond to the respective traits.}
  \item{abundances}{A matrix containing the abundances of all species for all plots. Rows correspond to species, columns to plots.}
}
\author{ Thijs Janzen & Fons van der Plas}

\examples{
Artificial.Data <- generate.Artificial.Data(n_species=40, n_traits=3, n_communities = 5,
                   occurence_distribution = 0.2, average_richness = 0.5, sd_richness = 0.2,
                   mechanism_random = FALSE)
}
