\name{data-treated.dataframe}
\alias{treated.dataframe}
\title{
Edges from the genetic interaction network created by Bandyopadhyay et al. (2010) using yeast treated with the DNA-damaging agent methyl methanesulfonate (MMS) 
}                                                                        
\description{ 
	It has been demonstrated that some genetic interactions (GIs) are condition-dependant. This data frame contains a filtered set of the interactions identified in yeast when exposed to DNA-damaging conditions. The data frame contains 415 genes linked by 4428 positive and negative GIs. GIs with a score between -1.75 and 1.75 have been removed.
	
	The data frame contains 3 columns, two containing the binding partners and one containing the GI scores. 
}
\references{
Bandyopadhyay, S., Mehta, M., Kuo, D., Sung, M. K., Chuang, R., Jaehnig, E. J., Bodenmiller, B., Licon, K., Copeland, W., Shales, M., Fiedler, D., Dutkowski, J., Guenole, A., Attikum, H., Shokat, K. M., Kolodner, R. D., Huh, W. K., Aebersold, R., Keogh, M. C., Krogan, N. and Ideker, T. (2010) \emph{Rewiring of genetic networks in response to DNA damage} Science, 330: 1385-1389. 
}
\examples{
# Load the data frame, create the network and display the number of vertices and edges 
# contained within.
data(treated.dataframe)
g.treated <- graph.data.frame(treated.dataframe, directed=FALSE)
vcount(g.treated)
ecount(g.treated)
}
