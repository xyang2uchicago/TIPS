library(igraph)
library(dplyr)
library(data.table)
library(qlcMatrix)
library(foreach)
library(doParallel)
library(MLmetrics)

########## BEGINNING OF USER INPUT ##########

wd = "/project/xyang2/felixy/GSE87038_weighted/" # Run on HPC
setwd(paste0(wd, "results/"))

db <- "GSE87038"

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

inputDir <- "PPI_weight/"
outputDir <- "PPI_weight/"

s <- "combined" # specificity method
failureAnalysis <- TRUE
V_attack <- TRUE # Targeted vertex attack
E_attack <- TRUE # Targeted edge attack

N_vertex <- 1000 # Number of targeted vertex attack simulations per network
N_edge <- 100 # Number of targeted vertex attack simulations per network

########## END OF USER INPUT ##########

### method 1 -- Random Protein failures in the same condition #######################################
# Think of it like: Testing how the same patient's cells respond to 100 different random genetic damages
# Goal: To assesses vulnerability to random failures
# Method: uisng robustness_MonteCarlo(..., measure='random')
# Network: Original graph (unchanged)
# Attack: Random vertex removal
# Question: "How robust is THIS network to random failures?"
# output: 100 simulations of random removal orders on the same network → averaged into single curve
# expect:  higher robustness, low variability
# interpretation: higher score -> genes are consistent resilient in this network
####################################################

# refer to 11.2.0_weighted_graph_attack_robustness.R
file <- paste0(inputDir, paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds"))
graph_list <- readRDS(file)
if (failureAnalysis) {
    cat('failureAnalysis in progress \n')
    tmp <- list()
    for (j in seq_along(graph_list)) {
        if (length(E(graph_list[[j]])) == 1) {
            tmp[[j]] <- NULL
        } else {
            tmp[[j]] <- robustness_MonteCarlo(graph_list[[j]], "edge", "random", N = 1e2)
        }
    }
    names(tmp) <- names(graph_list)
    failure.edge <- rbindlist(tmp, idcol = names(graph_list))
    saveRDS(failure.edge, file = paste0(outputDir, "failure.edge_100_simplified_", s, "weighted.rds"))


    failure.vertex <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", measure = "random", N = 1e2), idcol = names(graph_list))
    saveRDS(failure.vertex, file = paste0(outputDir, "failure.vertex_100_simplified_", s, "weighted.rds"))
    cat("failureAnalysis done \n")
}

vn <- pn <- array(dim = length(graph_list))
names(vn) <- names(pn) <- names(graph_list)
for (j in names(vn)) {
    vn[j] <- length(V(graph_list[[j]])) #  counts the number of vertices).
    pn[j] <- (graph_list[[j]] %>% degree() %>% mean()) / vn[j] # the average degree of the vertices in each graph, normalized by the number of vertices in the graph.
}

### method 2 Targeted Drug Attack on 100 Network Variants ###################################
# Think of it like: Testing how 100 different patients respond to the same precision cancer drug
# Goal: To explore how network structure affects robustness to strategic attacks
# Method: with rewiring + targeted attack strategy
# Network: Randomly rewired variants of original
# Attack: Targeted (highest betweenness first)
# Question: "How does robustness vary across similar network topologies?"
# output: 100 simulations of targeted attacks on 100 different network architectures → 100 distinct robustness measurements
# expect  lower robustness, higher variability
# interpretation: higher score -> genes are connected in a vulnerable structure
####################################################

if (V_attack) {
    cat('V_attack simulation in progress \n')
    N = N_vertex
    set.seed(1234)
    attac_V_random <- matrix(nrow = N, ncol = length(graph_list))
    colnames(attac_V_random) <- names(graph_list)
    for (j in colnames(attac_V_random)) {
        g <- graph_list[[j]]
        for (i in 1:N) {
            # This line of code rewires the graph g (from the graph_list) with a probability pn[j] for each edge,
            # where pn[j] is the average degree of the graph normalized by the number of vertices.
            # Node and Edge Count Remain Constant, and the node  with the largest betweenness after reviring is then removed
            # In contrast, robustness_MonteCarlo(g, 'vertex', 'random') randomly select node to remove
            g <- rewire(g, each_edge(prob = pn[j])) # rewires the graph edges based on the probability stored in pn[j]. It rewires the endpoints of the edges with a constant probability uniformly randomly to a new vertex in a graph.
            res <- robustness_MonteCarlo(g, "vertex", "btwn.cent") # calculate the betweenness centrality
            attac_V_random[i, j] <- Area_Under_Curve(res$removed.pct, res$comp.pct)
        }
    }
    saveRDS(attac_V_random, file = paste0(outputDir, "AUC_compt.pct_attac_V_random_", N, "runs_", s, "weighted.RDS"))
    cat("V_attack done \n")
}

if (E_attack) {
    cat('E_attack simulation in progress \n')
    N = N_edge
    # To save the results
    attac_E_random <- matrix(nrow = N, ncol = length(graph_list))
    colnames(attac_E_random) <- names(graph_list)

    set.seed(1234)

    for (j in colnames(attac_E_random)) {
        g <- graph_list[[j]]
        if (length(E(g)) > 1) {
            for (i in 1:N) {
                g <- rewire(g, igraph::each_edge(prob = pn[j]))

                if ("weight" %in% igraph::edge_attr_names(g)) {
                    # igraph::edge_betweenness() uses distance graph weights, but E(g) uses connection weights, thus we invert it.
                    E(g)$btwn <- edge_betweenness(g, weights = 1/E(g)$weight)
                } else {
                    E(g)$btwn <- edge_betweenness(g) # otherwise encountering the NA value or the original betweenness was retained
                }
                res <- robustness_MonteCarlo(g, "edge", "btwn.cent", N)
                attac_E_random[i, j] <- Area_Under_Curve(res$removed.pct, res$comp.pct)
            }
        }
    }

    saveRDS(attac_E_random, file = paste0(outputDir, "AUC_compt.pct_attac_E_random_", N, "runs_", s, "weighted.RDS"))
    cat("E_attack done \n")
}
