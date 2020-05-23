
rm(list = ls())

# loads required packages 
for (pkg in c("tidyverse", "igraph", "backbone", "reshape2")) {library(pkg, character.only = TRUE)}

# loading data 
setwd("/sfs/qumulo/qhome/kb7hp/data/epik-data")
raw_ev_edgelist <- read_csv("052918 - Important Hearings - Unimodal - Cuts Across Years - Evidence.csv")

raw_ev_edgelist
# We can see this is just our original edgelist.
# Our original ergms did not use weights, as statnet defaults to a network with no weights. 
# To get weights, let's run a group_by on this edgelist. 

ev_edgelist <- raw_ev_edgelist %>% 
  rename(from = Source, to = Target) %>% 
  group_by(from, to) %>% # collapses all duplicate rows 
  summarise(weight = n()) %>% # gives the weight for duplicate rows 
  arrange(-weight) %>% # orders by edge weight 
  ungroup(); ev_edgelist

# And now we see this weighted edgelist and we see some issues: 
# First, you see that the relations are duplicated (WilliamHDietz-MargaretHamburg and MargaretHamburg-WilliamHDietz) both show up. 
# Second, we see that the edge weights are in the thousands, which means that these actors are sharing evidence 2000+ times 
# over the course of the 786 policy documents we examined. Is that right? It doesn't seem plausible to me. 
# In our non-weighted, neiether of these issues are a concern since we only factored in whether evidence was shared or not. 
# Regardless, addressing the first issue (i.e. duplicate rows) can be fixed by using igraph to collapse the relations. 

ev_network <- simplify(graph.data.frame(ev_edgelist, directed = FALSE), 
                        remove.loops = FALSE, 
                        edge.attr.comb = igraph_opt("edge.attr.comb"))
is_weighted(ev_network)

final_edgelist <- reshape2::melt(as.matrix(as_adjacency_matrix(ev_network, type = "lower", attr = "weight", sparse = T)))

final_edgelist <- final_edgelist %>% 
  rename(Source = Var1, Target = Var2) %>% 
  mutate(Weight = round(value / 2,0)) %>% 
  filter(value != 0) %>% 
  select(Source, Target, Weight) %>% 
  arrange(-Weight); final_edgelist

# Here, you can see that the final product is a weighted edgelist where the duplicates are removed 
# If Matt feels that these weights are within the realm of possibility than we can start running ERGMs
# If these weights are not within the realm of possibility, 
# we need to go back to the original data output to find out why these are not plausible. 

setwd("/sfs/qumulo/qhome/kb7hp/data/epik-data")
write_csv(final_edgelist, "EPIK_Evidence_WeightedEdgelist_052320.csv")

