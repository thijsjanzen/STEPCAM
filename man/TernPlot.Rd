\name{TernPlot}
\alias{TernPlot}
\title{
  Create a ternary plot from results of the STEPCAM_ABC function
}
\description{
  This function plots the obtained combinations of dispersal assembly, habitat filtering and limiting similarity of best fitting models in a Ternary plot, similar to Figure 3 in van der Plas et al. 2015
}
\usage{
TernPlot(output)
}
\arguments{
  \item{output}{A list with the following elements: \cr
           - \code{DA}  a vector containing the posterior distribution for the parameter dispersal assembly \cr
           - \code{HF}  a vector containing the posterior distribution for the parameter habitat filtering \cr
           - \code{LS}  a vector containing the posterior distribution for the parameter limiting similarity \cr
            This list is generated as output by the function STEPCAM_ABC}
}

\references{Fons van der Plas, Thijs Janzen, Alejandro Ordonez, Wimke Fokkema, Josephine Reinders, Rampal S. Etienne, and Han Olff 2015. A new modeling approach estimates the relative importance of different community assembly processes. Ecology 96:1502-1515.}
\author{Fons van der Plas. Adapted from the ternaryplot function from the 'vcd' package}
\examples{
Artificial.Data <- generate.Artificial.Data(numSpecies=40,numTraits=3,numCommunities=5,
                   occurence_distribution=0.2,average_richness=0.5,SD_richness=0.2,
                   random.Mechanism=FALSE);
O <- STEPCAM_ABC(Artificial.Data$abundances,Artificial.Data$traits,
                 numParticles=10,n_traits=3,plot_number=1,stopRate=0.8);
TernPlot(O);
}
