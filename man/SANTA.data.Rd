\name{SANTA.data}
\alias{edgelist.humannet}
\alias{edgelist.intact}
\alias{g.bandyopadhyay.treated}
\alias{g.bandyopadhyay.untreated}
\alias{g.costanzo.cor}
\alias{g.costanzo.raw}
\alias{g.srivas.high}
\alias{g.srivas.untreated}
\alias{go.entrez}
\alias{rnai.cheung}
\title{
Pre-processed dataset for the SANTA vignette
}
\description{
Pre-processed network and expression data for use in the accompanying vignette. 
}
\details{
\describe{
  \item{edgelist.humannet}{Data frame of functional interactions from the HumanNet network with edge distances derived from the accompanying log-likelihood scores.}
  \item{edgelist.intact}{Data frame of \emph{H. sapiens} physical interactions from the IntAct database downloaded on 2013-05-02.}
  \item{g.bandyopadhyay.treated}{\emph{S. cerevisiae} interaction network created from the MMS-treated GI data from the study by Bandyopadhyay et al. Interactions represent correlation in GI profile. Network is an \code{igraph} object.}
  \item{g.bandyopadhyay.untreated}{\emph{S. cerevisiae} interaction network created from the untreated GI data from the study by Bandyopadhyay et al. Interactions represent correlation in GI profile. Network is an \code{igraph} object.}
  \item{g.costanzo.cor}{\emph{S. cerevisiae} interaction network created from the GI data from the study by Costanzo et al. Interactions represent correlation in GI profile. Network is an \code{igraph} object.}
  \item{g.costanzo.raw}{\emph{S. cerevisiae} interaction network created from the GI data from the study by Costanzo et al. Interactions represent raw GIs. Network is an \code{igraph} object.}
  \item{g.srivas.high}{\emph{S. cerevisiae} interaction network created from the high UV-dosage GI data from the study by Srivas et al. Interactions represent raw GIs. Network is an \code{igraph} object.}
  \item{g.srivas.untreated}{\emph{S. cerevisiae} interaction network created from the untreated GI data from the study by Srivas et al. Interactions represent raw GIs. Network is an \code{igraph} object.}
  \item{go.entrez}{Genes associated with the GO term GO:0042981 (regulation of apoptotic process).}
  \item{rnai.cheung}{Matrix of gene-wise essentiality p-values for 6 cancer cell lines. Rows represent genes and columns represent cancers. This data was created from the RNAi screens conducted by Cheung et al. The weight of evidence approach was used to compute essentiality scores for each shRNA and GENE-E was used to collapse the shRNA scores into gene-wise scores.}
}
}
\references{
Lee, I., Blom, U.M., Wang, P.I. et al. (2011). \emph{Prioritizing candidate disease genes by network-based boosting of genome-wide association data.} Genome Research. 21, 1109-21. 

Orchard, S., Ammari, M., Aranda, B. et al. (2014). \emph{The MIntAct project - IntAct as a common curation platform for 11 molecular interaction databases.} Nucleic Acids Research. 42:1, D358-63.

Costanzo, M., Baryshnikova., A., Bellay, J. et al. (2010). \emph{The Genetic Landscape of a Cell.} Science. 327:5964, 425-31. 

Bandyopadhyay, S., Mehta, M., Kuo, D. et al. (2010) \emph{Rewiring of Genetic Networks in Response to DNA Damage.} Science. 330, 1385-90.

Srivas, R., Costelloe, T., Carvunis, A. et al. (2013) \emph{A UV-induced genetic network links the RSC complex to nucleotide excision repair and shows dose-dependant rewiring.} Cell Reports. 5:6, 1714-24.

Ashburner, M., Ball, C.A., Blake, J.A. et al. (2000) \emph{Gene Ontology: tool for the unification of biology.} Nature Genetics. 25, 25-9.

Cheung, H.W, Cowley, G.S., Weir, B.A. et al. (2011) \emph{Systematic investigation of genetic vulnerabilities across cancer cell lines reveals lineage-specific dependencies in ovarian cancer.} PNAS. 108:30, 12372-7.
}
\examples{
data(g.bandyopadhyay.treated)
}
