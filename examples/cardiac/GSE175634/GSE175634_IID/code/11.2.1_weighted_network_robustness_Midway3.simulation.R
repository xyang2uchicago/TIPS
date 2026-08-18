## by setting different working direction -- setwd(<...> , we are testing the impact of <shrink=TRUE> or <shrink=FALSE>

# -- simialr to old 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R 
# -- but calling robustness_MonteCarlo() which is based on brainGraph::robustness function
# This function performs a “targeted attack” of a graph or a “random failure” analysis, calculating the size of the largest component after edge or vertex removal.

# In a targeted attack, it will sort the vertices by either degree or betweenness centrality (or sort edges by betweenness), and successively remove the top vertices/edges. Then it calculates the size of the largest component.

# In a random failure analysis, vertices/edges are removed in a random order.
## we added a new function new_centr_betw () called by robustness_MonteCarlo() in the source code file 'celltype_specific_weight_v2.R'
#
# robustness_MonteCarlo() is an revised version of brainGraph::robustness()
# "This function analyzes network robustness by systematically removing vertices or edges 
# and measuring how the size of the largest connected component changes over time. 
# It supports both targeted attacks (removing high centrality nodes/edges first) and 
# random failure simulations (Monte Carlo approach with N repetitions). 
# The output shows the network's vulnerability profile, indicating how quickly the 
# network fragments under different attack strategies."


# rcchelp balance  
# rcchelp usage
# squeue -u xyang2
# squeue -p bigmem --state=PD | wc -l
# squeue -p caslake --state=PD | wc -l
# squeue -p gpu --state=PD | wc -l
# squeue -p amd-hm --state=PD | wc -l

# sinteractive -p gpu --account=pi-xyang2 --gres=gpu:1 --mem=180GB  --time=8:00:00 -c 4 

# scp -p -r  F:/projects/scRNA/source/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9/11.2.1_weighted_network_robustness_Midway3.simulation.R   xyang2@midway3.rcc.uchicago.edu:/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/.
source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v10.R')
# module load python/anaconda-2022.05
# source activate /project/xyang2/software-packages/env/velocity_2025Feb_xy  
   
library(igraph)
library(dplyr)
library(data.table)
library(brainGraph)

library(foreach)
library(doParallel)
# n_cores = 6

# setwd('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weight/')
setwd('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_IID/')
#source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/robustness2.R')

source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP_update_06162025.R')
source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v10.R')

### method 1 -- Random Protein failures in the same condition #######################################
# Think of it like: Testing how the same patient's cells respond to 100 different random genetic damages
# Goal: To assesses vulnerability to random failures
# Method: uisng robustness_MonteCarlo(..., measure='random')
# Network: Original graph (unchanged)
# Attack: Random vertex removal
# Question: "How robust is THIS network to random failures?"
# output: 100 simulations of random removal orders on the same network → averaged into single curve
# expect:  higher robustness, low variability
# interpretation: higher score -> genes are consistent reilient in this network
####################################################

# refer to 11.2.0_weighted_graph_attack_robustness.R
 #for(s in c("combined","ratio", "zscore", "diff")){
s = "combined"
file = paste0('GSE175634_IID_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS( file)  

failureAnalysis=TRUE
if(failureAnalysis){
	tmp = list()
	for(j in 1:length(graph_list)){
		if(length(E(graph_list[[j]]))==1)  tmp[[j]] =NULL else
		  tmp[[j]] = robustness_MonteCarlo(graph_list[[j]],'edge', 'random', N=1e2)
	} 
	names(tmp) = names(graph_list)
	failure.edge <- rbindlist(tmp, idcol=names(graph_list))
	saveRDS(failure.edge, file= paste0('failure.edge_100_simplified_',s,'weighted.rds'))


	failure.vertex <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", measure='random', N=1e2), idcol=names(graph_list))
	saveRDS(failure.vertex, file= paste0('failure.vertex_100_simplified_',s,'weighted.rds'))
}
   	#### calcualte empritical p-valeus on midway3
	## refer to 

library(MLmetrics)
library(sm)

vn = pn = array(dim=length(graph_list))
names(vn) = names(pn) = names(graph_list)
for(j in names(vn)){
	vn[j] = length(V(graph_list[[j]]))  #  counts the number of vertices).
	pn[j] = (graph_list[[j]] %>% degree %>% mean)/vn[j] # the average degree of the vertices in each graph, normalized by the number of vertices in the graph.
}
cat('failureAnalysis done \n')

N = 1000

### method 2 Targeted Drug Attack on 100 Network Variants ###################################
# Think of it like: Testing how 100 different patients respond to the same precision cancer drug
# Goal: To explore how network structure affects robustness to strategic attacks
# Method: with rewiring + targeted attack strategy
# Network: Randomly rewired variants of original
# Attack: Targeted (highest betweenness first)
# Question: "How does robustness vary across similar network topologies?"
# output: 100 simulations of targeted attacks on 100 different network architectures → 100 distinct robustness measurements
# expect  lower robustness, higher variability
# interpretation: higher score -> genes are connected in a vulenrable structure
####################################################

V_attack=TRUE
if(V_attack){
	set.seed(1234)
	attac_V_random = matrix(nrow=N, ncol=length(graph_list))  
    colnames(attac_V_random) = names(graph_list)
	for(j in colnames(attac_V_random)){
		g0 = graph_list[[j]]
		for(i in 1:N){   
			# This line of code rewires the graph g (from the graph_list) with a probability pn[j] for each edge, 
			# where pn[j] is the average degree of the graph normalized by the number of vertices.
			# Node and Edge Count Remain Constant, and the node  with the largest betweenness after reviring is then removed
			# In contrast, robustness_MonteCarlo(g, 'vertex', 'random') randomly select node to remove
		  g = rewire(g0, each_edge(prob = pn[j]))  # rewires the graph edges based on the probability stored in pn[j]. It rewires the endpoints of the edges with a constant probability uniformly randomly to a new vertex in a graph.
		  res = robustness_MonteCarlo(g, 'vertex', 'btwn.cent') # calculate the betweenness centrality 
		  attac_V_random[i,j] = Area_Under_Curve(res$removed.pct, res$comp.pct)
		 }
	}
	saveRDS(attac_V_random, file = paste0('AUC_compt.pct_attac_V_random_',N,'runs_',s,'weighted.RDS'))
}
cat('V_attack done \n')
	
E_attack=TRUE
if(E_attack) {
    N <- 100
    # To save the results
    attac_E_random <- matrix(nrow=N, ncol=length(graph_list))   
    colnames(attac_E_random) <- names(graph_list)
  
    set.seed(1234)  

	for(j in colnames(attac_E_random)){
	# j = 16 
		g = graph_list[[j]]
		if( length(E(g))>1) for(i in 1:N){      
		  g = rewire(g, igraph::each_edge(prob = pn[j]))
		  
		  if ("weight" %in% igraph::edge_attr_names(g)) {
		    E(g)$btwn = edge_betweenness(g, weights = 1/E(g)$weight)
		  } else {
			E(g)$btwn = edge_betweenness(g)  # otherwise encountering the NA value or the original betweenness was retained
		  }
		 res = robustness_MonteCarlo(g, 'edge', 'btwn.cent', N)
		 attac_E_random[i,j] = Area_Under_Curve(res$removed.pct, res$comp.pct)
		 }
	}

	  
  saveRDS(attac_E_random, file = paste0('AUC_compt.pct_attac_E_random_',N,'runs_',s,'weighted.RDS'))
}
cat('E_attack done \n')

# scp -p -r  xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_IID/*  F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID/.

  