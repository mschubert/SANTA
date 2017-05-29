## ----setup, echo=FALSE---------------------------------------------------
options(width = 75)
options(useFancyQuotes=FALSE)


## ----hook-printfun, echo=FALSE-------------------------------------------
library(knitr)
library(formatR)
knit_hooks$set(printfun = function(before, options, envir) {
    if (before) return()
    txt = capture.output(dump(options$printfun, '', envir = envir))
    ## reformat if tidy=TRUE
    if (options$tidy)
        txt = tidy.source(text=txt, output=FALSE,
                          width.cutoff=30L, keep.comment=TRUE,
                          keep.blank.line=FALSE)$text.tidy
    paste(c('\n```r\n', txt, '\n```\n'), collapse="\n")
    })

## ----case_i_graph_create-------------------------------------------------
# generate the simulated network
require(SANTA)
require(igraph)
set.seed(1) # for reproducibility
g <- barabasi.game(n=500, power=1, m=1, directed=F) 

## ----case_i_compute_d----------------------------------------------------
# measure the distance between pairs of vertices in g
dist.method <- "shortest.paths"
D <- DistGraph(g, dist.method=dist.method, verbose=F) 

## ----case_i_compute_b----------------------------------------------------
# place the distances into discreet bins
B <- BinGraph(D, verbose=F)

## ----case_i_trials-------------------------------------------------------
cluster.size <- 5
s.to.use <- c(10, 20, 50, 100, 500)
n.trials <- 10
pvalues <- array(0, dim=c(n.trials, length(s.to.use)), 
    dimnames=list(NULL, as.character(s.to.use)))

# run the trials for each value of s
for (s in s.to.use) {    	
    for (i in 1:n.trials) {
        # generate the hit set
        seed.vertex <- sample(vcount(g), 1) # identify seed
        sample.set <- order(D[seed.vertex, ])[1:s]
        hit.set <- sample(sample.set, cluster.size)
        
        # measure the stength of association
        g <- set.vertex.attribute(g, name="hits", 
            value=as.numeric(1:vcount(g) %in% hit.set))
        pvalues[i, as.character(s)] <- Knet(g, nperm=100, 
            dist.method=dist.method, vertex.attr="hits", 
            B=B, verbose=F)$pval
    }
}

## ----case_i_plot, results="asis", echo=FALSE-----------------------------
boxplot(-log10(pvalues), xlab="cutoff", ylab="-log10(p-value)")

## ----case_i_cleanup, echo=FALSE------------------------------------------
# cleanup
rm(B, D, g)

## ----case_ii_create_network----------------------------------------------
# create the network
n.nodes <- 12
edges <- c(1,2, 1,3, 1,4, 2,3, 2,4, 3,4, 1,5, 5,6, 
           2,7, 7,8, 4,9, 9,10, 3,11, 11,12)
weights1 <- weights2 <- rep(0, n.nodes)
weights1[c(1,2)] <- 1
weights2[c(5,6)] <- 1 

g <- graph.empty(n.nodes, directed=F)
g <- add.edges(g, edges)
g <- set.vertex.attribute(g, "weights1", value=weights1)
g <- set.vertex.attribute(g, "weights2", value=weights2)

## ----case_ii_create_plot, results="asis", echo=FALSE---------------------
par(mfrow=c(1,2))

colors <- rep("grey", n.nodes)
colors[which(weights1 == 1)] <- "red"
g <- set.vertex.attribute(g, "color", value=colors)
plot(g)

colors <- rep("grey", n.nodes)
colors[which(weights2 == 1)] <- "red"
g <- set.vertex.attribute(g, "color", value=colors)
plot(g)

par(mfrow=c(1,1))
g <- remove.vertex.attribute(g, "color")

## ----case_ii_pvalues-----------------------------------------------------
# set 1
Knet(g, nperm=100, vertex.attr="weights1", verbose=F)$pval
Compactness(g, nperm=100, vertex.attr="weights1", verbose=F)$pval

# set 2
Knet(g, nperm=100, vertex.attr="weights2", verbose=F)$pval
Compactness(g, nperm=100, vertex.attr="weights2", verbose=F)$pval

## ----case_iii_load_networks----------------------------------------------
# load igraph objects
data(g.costanzo.raw)
data(g.costanzo.cor)
networks <- list(raw=g.costanzo.raw, cor=g.costanzo.cor)
network.names <- names(networks)
network.genes <- V(networks$raw)$name 
# genes identical across networks

## ----case_iii_go---------------------------------------------------------
# obtain the GO term associations from org.Sc.sgd.db package
library(org.Sc.sgd.db)
xx <- as.list(org.Sc.sgdGO2ALLORFS)
go.terms <- c("GO:0000082", "GO:0003682", "GO:0007265", 
              "GO:0040008", "GO:0090329")
# apply the GO terms to the networks
for (name in network.names) {
    for (go.term in go.terms) {
        networks[[name]] <- set.vertex.attribute(
            networks[[name]], name=go.term, 
            value=as.numeric(network.genes %in% xx[[go.term]]))
    }
} 

## ----case_iii_association------------------------------------------------
# results <- list()
# for (name in network.names) {
#     results[[name]] <- Knet(networks[[name]], nperm=1000, 
#         vertex.attr=go.terms, edge.attr="distance", verbose=F)
#     results[[name]] <- sapply(results[[name]], 
#         function(res) res$pval)
# }

## ----case_iii_plot-------------------------------------------------------
# p.values <- array(unlist(results), dim=c(length(go.terms), 
#       length(network.names)), dimnames=list(go.terms, 
#       network.names))
# p.values.ml10 <- -log10(p.values)
# axis.range <- c(0, max(p.values.ml10)) 
# plot(p.values.ml10[, "raw"], p.values.ml10[, "cor"], asp=1, 
#       xlim=axis.range, ylim=axis.range, bty="l", 
#       xlab="-log10 of the p-value in the raw GI network", 
#       ylab="-log10 of the p-value in the correlation network", 
#       main="")
# abline(0, 1, col="red")

## ----case_iii_cleanup, echo=FALSE----------------------------------------
# cleanup
rm(g.costanzo.cor, g.costanzo.raw, network.genes, networks, xx)

## ----case_iv_load_igraph-------------------------------------------------
# load igraph objects
data(g.bandyopadhyay.treated)
data(g.bandyopadhyay.untreated)
networks <- list(
    treated=g.bandyopadhyay.treated, 
    untreated=g.bandyopadhyay.untreated
)
network.names <- names(networks)

## ----case_iv_go----------------------------------------------------------
# obtain GO term associations
library(org.Sc.sgd.db)
xx <- as.list(org.Sc.sgdGO2ALLORFS)
# change to use alternative GO terms
associated.genes <- xx[["GO:0006974"]] 
associations <- sapply(networks, function(g) 
    as.numeric(V(g)$name %in% associated.genes), 
    simplify=F)
networks <- sapply(network.names, function(name) 
    set.vertex.attribute(networks[[name]], "rdds", 
    value=associations[[name]]), simplify=F)

## ----case_iv_association-------------------------------------------------
# results <- sapply(networks, function(g) Knet(g, nperm=1000, 
#     dist.method="shortest.paths", vertex.attr="rdds", 
#     edge.attr="distance"), simplify=F)
# p.values <- sapply(results, function(res) res$pval)
# p.values

## ----case_iv_plot--------------------------------------------------------
# plot(results$treated)
# plot(results$untreated)

## ----case_iv_cleanup, echo=FALSE-----------------------------------------
# cleanup
rm(xx, networks, g.bandyopadhyay.treated, g.bandyopadhyay.untreated, associated.genes)

## ----case_v_load---------------------------------------------------------
# laod igraph object
data(g.srivas.high)
data(g.srivas.untreated)
networks <- list(
    high=g.srivas.high, 
    untreated=g.srivas.untreated
)
network.names <- names(networks)

## ----case_v_go-----------------------------------------------------------
# obtain GO term associations
library(org.Sc.sgd.db)
xx <- as.list(org.Sc.sgdGO2ALLORFS)
associated.genes <- xx[["GO:0000725"]] 
associations <- sapply(networks, function(g) 
    as.numeric(V(g)$name %in% associated.genes), 
    simplify=F)
networks <- sapply(network.names, function(name)
    set.vertex.attribute(networks[[name]], "dsbr", 
    value=associations[[name]]), simplify=F)

## ----case_v_association--------------------------------------------------
# p.values <- sapply(networks, function(g) 
#       Knet(g, nperm=1000, dist.method="shortest.paths", 
#       vertex.attr="dsbr", edge.attr="distance")$pval)
# p.values

## ----case_v_cleanup, echo=FALSE------------------------------------------
# cleanup
rm(xx, networks, g.srivas.high, g.srivas.untreated)

## ----case_vi_rnai--------------------------------------------------------
# import and convert RNAi data
data(rnai.cheung)
rnai.cheung <- -log10(rnai.cheung)
rnai.cheung[!is.finite(rnai.cheung)] <- max(rnai.cheung[is.finite(rnai.cheung)])
rnai.cheung[rnai.cheung < 0] <- 0

## ----case_vi_networks----------------------------------------------------
# import and create IntAct network
data(edgelist.intact)
g.intact <- graph.edgelist(as.matrix(edgelist.intact), 
    directed=FALSE)

# import data and create HumanNet network
data(edgelist.humannet)
g.humannet <- graph.edgelist(as.matrix(edgelist.humannet)[,1:2], 
    directed=FALSE)
g.humannet <- set.edge.attribute(g.humannet, "distance", 
    value=edgelist.humannet$distance)
networks <- list(intact=g.intact, humannet=g.humannet)

## ----case_vi_weights-----------------------------------------------------
network.names <- names(networks)
network.genes <- sapply(networks, get.vertex.attribute, 
    name="name", simplify=F)
rnai.cheung.genes <- rownames(rnai.cheung)
cancers <- colnames(rnai.cheung)

for (cancer in cancers) {
    for (name in network.names) {
        vertex.weights <-rep(NA, vcount(networks[[name]]))
        names(vertex.weights) <- network.genes[[name]]
        common.genes <- rnai.cheung.genes[rnai.cheung.genes
            %in% network.genes[[name]]] 
        vertex.weights[common.genes] <- rnai.cheung[common.genes, cancer]
        networks[[name]] <- set.vertex.attribute(networks[[name]], 
            cancer, value=vertex.weights)
    }
}

## ----case_vi_association-------------------------------------------------
#knet.res <- sapply(networks, Knet, nperm=100, 
#   dist.method="shortest.paths", vertex.attr=cancers, 
#   edge.attr="distance", simplify=F)
#p.values <- sapply(knet.res, function(i) sapply(i, 
#   function(j) j$pval))

## ----case_vi_cleanup, echo=FALSE-----------------------------------------
# cleanup
rm(rnai.cheung, rnai.cheung.genes, networks, network.genes, edgelist.humannet, edgelist.intact, common.genes, g.humannet, g.intact)

## ----case_vii_bionet-----------------------------------------------------
#library(BioNet)

## ----case_vii_network----------------------------------------------------
# # required parameters
# n.nodes <- 1000
# n.hits <- 10
# clusters <- 3
# 
# # create network and spread hits across 3 clusters
# g <- barabasi.game(n=n.nodes, power=1, m=1, directed=FALSE) 
# g <- SpreadHits(g, h=n.hits, clusters=clusters, distance.cutoff=12, 
#     lambda=10, dist.method="shortest.paths", verbose=FALSE)

## ----case_vii_pvalues----------------------------------------------------
# # simulate p-values
# library(msm)
# hits <- which(get.vertex.attribute(g, "hits") == 1)
# p.values <- runif(vcount(g))
# names(p.values) <- as.character(1:vcount(g))
# p.values[as.character(hits)] <- rtnorm(n.hits * clusters, mean=0, 
#     sd=10e-6, lower=0, upper=1)

## ----case_vii_apply_bionet-----------------------------------------------
# # apply BioNet to the network and p-values
# bum <- fitBumModel(p.values, plot=F)
# scores <- scoreNodes(network=g, fb=bum, fdr=0.1)
# module <- runFastHeinz(g, scores)

## ----case_vii_apply_knode------------------------------------------------
# # apply Knode to the network
# g <- set.vertex.attribute(g, name="pheno", value=-log10(p.values))  
# knode.results <- Knode(g, dist.method="diffusion", 
#     vertex.attr="pheno", verbose=FALSE)

## ----case_vii_numbers----------------------------------------------------
# # number of hits identified by BioNet 
# sum(hits %in% as.numeric(V(module)$name))
# 
# # number of hits ranked within the top 30 by Knode
# sum(hits %in% as.numeric(names(knode.results)[1:(n.hits * clusters)]))

## ----case_vii_plot-------------------------------------------------------
# # create subnetworks
# g.bionet <- g.knode <- induced.subgraph(g, hits)
# color.bionet <- color.knode <- rep("grey", vcount(g.bionet))
# color.bionet[hits %in% as.numeric(V(module)$name)] <- "blue"
# color.knode[hits %in% as.numeric(names(knode.results)[1:(n.hits * clusters)])] <- "green"
# g.bionet <- set.vertex.attribute(g.bionet, "color", value=color.bionet)
# g.knode <- set.vertex.attribute(g.knode, "color", value=color.knode)
# 
# # plot subnetworks
# par(mfrow=c(1,2))
# plot(g.bionet)
# plot(g.knode)

## ----case_vii_cleanup, echo=FALSE----------------------------------------
# # cleanup
# rm(vertex.weights, p.values, g)

## ----case_viii_load------------------------------------------------------
# # load required package
# library(SANTA)
# library(BioNet)
# library(DLBCL)
# data(exprLym) 
# data(dataLym) 
# data(interactome) 

## ----case_viii_pvalues---------------------------------------------------
# # extract entrez IDs
# library(stringr)
# 
# # aggregate p-values
# pvals <- cbind(dataLym$t.pval, dataLym$s.pval)
# pval <- aggrPvals(pvals, order=2, plot=F)
# names(pval) <- dataLym$label

## ----case_viii_network---------------------------------------------------
# # derive Lymphochip-specific network
# network <- subNetwork(featureNames(exprLym), interactome)
# network <- largestComp(network) # use only the largest component
# network <- igraph.from.graphNEL(network, name=T, weight=T) 
# network <- simplify(network)

## ----case_viii_run_bionet------------------------------------------------
# # run BioNet on the Lymphochip-specific network and aggregate p-values
# fb <- fitBumModel(pval, plot=F)
# scores <- scoreNodes(network, fb, fdr=0.001)
# module <- runFastHeinz(network, scores)
# extract.entrez <- function(x) str_extract(str_extract(x, 
#   "[(][0-9]+"), "[0-9]+")
# bionet.genes <- extract.entrez(V(module)$name)

## ----case_viii_run_knode-------------------------------------------------
# # convert p-values to vertex weights
# vertex.weights <- -log10(pval)[get.vertex.attribute(network, 
#   "name")]
# network <- set.vertex.attribute(network, name="pheno", 
#   value=vertex.weights)
# 
# # run Knode on the Lymphochip-specific network and 
# # converted aggregate p-values
# knode.results <- Knode(network, dist.method="diffusion", 
#     vertex.attr="pheno", verbose=F)
# knode.genes <- extract.entrez(names(knode.results)[1:vcount(module)])

## ----case_viii_compare---------------------------------------------------
# data(go.entrez)
# sum(go.entrez %in% bionet.genes)
# sum(go.entrez %in% knode.genes)

## ----case_viii_cleanup, echo=FALSE---------------------------------------
# # cleanup
# rm(dataLym, exprLym, interactome, network, pval, pvals)

## ----session_info--------------------------------------------------------
sessionInfo()

## ----setup_references, echo=FALSE, results="hide"------------------------
suppressMessages(library("knitcitations"))
cleanbib()
suppressMessages(bib <- read.bibtex("refs.bib"))
citations <- c("White2003", "Costanzo2010", "Gaetan2010", "GeneOntology2000", "Beisser2010", "Bandyopadhyay2010", "Peri2003", "Cheung2011", "Orchard2014", "Lee2011", "Srivas2013a")
suppressWarnings(suppressMessages(citep(bib[citations])))

## ----print_references, results='asis', echo=FALSE------------------------
bibliography()

