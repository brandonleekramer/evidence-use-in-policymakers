
#########################################################################
# Section 1: Housekeeping and network construction    
#########################################################################

# loading packages 
rm(list = ls())

for (pkg in c("tidyverse", "statnet", "ggplot2", "ggnetwork",
              "ergm", "btergm", "statnet", "parallel")) {library(pkg, character.only = TRUE)}

setwd("/sfs/qumulo/qhome/kb7hp/data/epik-data")
EVDATAN <- read_csv("EPIK_Evidence_WeightedEdgelist_052320.csv")
EVDATAN <- EVDATAN %>% filter(Source != Target)

# Constructing the network 
ev_edgelist <- as.matrix(as.data.frame(EVDATAN %>% select(-Weight))) 
EVNETN=as.network(ev_edgelist, matrix.type="edgelist",directed=FALSE)
EVDEGN=degree(EVNETN)
EVNETN%e%"weight" <- EVDATAN$Weight 
plot.network(EVNETN, edge.lwd = get.edge.value(EVNETN,"weight"))
summary.network(EVNETN, print.adj = FALSE)

# Loading network attributes 
attrEVNETN <- read.csv(file.choose(), header = TRUE) 
EVNETN%v%"affil"<-attrEVNETN$AffiliationCategory # Affiliation category (1 = House; 2 = Senate; 3 = NGOs; 4 = Education Advocacy Organizations; Other); 5 = Private Organizations; 6 = New Media Organizations; 7 = Health Advocacy Organizations; 8 = Federal Government, Other; 9 = Professional Organizations and Lobbying Organizations; Missing = 999)
EVNETN%v%"gender"<-attrEVNETN$Gender # Gender of policymaker 
EVNETN%v%"political"<-attrEVNETN$PoliticalParty # Political party (when a governmental/political actor)
EVNETN%v%"polnotpol" <-attrEVNETN$Polnotpol # Proxy for having a governmental/political position or not 
EVNETN%v%"tenure"<-attrEVNETN$TenureTotal # Tenure in political position 
EVNETN%v%"branch" <-attrEVNETN$Branch # House or Senate 
# Five most prominent governmental sub-committees in shaping childhood obesity legislation from 2000-2014 
EVNETN%v%"appropriations" <-attrEVNETN$SenAppComm 
EVNETN%v%"help" <-attrEVNETN$SenHELPComm
EVNETN%v%"agriculture" <-attrEVNETN$SenAgComm
EVNETN%v%"energycomm" <-attrEVNETN$HouseEnergyComm
EVNETN%v%"educ_wf" <-attrEVNETN$HouseEducWFComm 

# Checking the number of cores available in the container 
detectCores() # Brandon drew from 16 cores in final models, but developed the models using an instance with 32 cores 
rivanna_cores <- 6

# LINKS 

# https://docs.google.com/spreadsheets/d/1UQ3cAHHzK6_P-0Rg06c3pvj0oyIbcqr0KZB6Qz4lL7E/edit#gid=1089563747


# BASELINE MODEL                                                                                                         BASELINE

ev <- ergm(EVNETN~edges, control=control.ergm(parallel=rivanna_cores, parallel.type="PSOCK"))
summary(ev)

# HOMOPHILY                                                                                                              HOMOPHILY
ev_aff <- ergm(EVNETN~edges + nodematch("affil"), 
                 control=control.ergm(rivanna_cores, 
                 MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_aff)
mcmc.diagnostics(ev_aff)
ev_affil_gof<-gof(ev_aff)
plot(ev_aff)

#gender 
ev_gndr <- ergm(EVNETN~edges + nodematch("gender"), 
                 control=control.ergm(rivanna_cores, 
                 MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gndr)
mcmc.diagnostics(ev_gndr)
ev_gndr_gof<-gof(ev_gndr)
plot(ev_gndr)

#political
ev_pol <- ergm(EVNETN~edges + nodematch("political"), 
                 control=control.ergm(rivanna_cores, 
               MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_pol)
mcmc.diagnostics(ev_pol)
ev_pol_gof<-gof(ev_pol)
plot(ev_pol)

#tenure
ev_ten <- ergm(EVNETN~edges + nodematch("tenure"), 
                 control=control.ergm(rivanna_cores, 
               MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_ten)
mcmc.diagnostics(ev_ten)
ev_ten_gof<-gof(ev_ten)
plot(ev_ten)

#branch
ev_branch <- ergm(EVNETN~edges + nodematch("branch"), 
               control=control.ergm(rivanna_cores, 
               MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_branch)
mcmc.diagnostics(ev_branch)
ev_branch_gof<-gof(ev_branch)
plot(ev_branch)

#appropriations
ev_app <- ergm(EVNETN~edges + nodematch("appropriations"), 
               control=control.ergm(rivanna_cores, 
               MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_app)
mcmc.diagnostics(ev_app)
ev_app_gof<-gof(ev_app)
plot(ev_app)

#help
ev_help <- ergm(EVNETN~edges + nodematch("help"), 
               control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_help)
mcmc.diagnostics(ev_help)
ev_help_gof<-gof(ev_help)
plot(ev_help)

#agriculture
ev_ag <- ergm(EVNETN~edges + nodematch("agriculture"), 
               control=control.ergm(rivanna_cores,  MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_ag)
mcmc.diagnostics(ev_ag)
ev_ag_gof<-gof(ev_ag)
plot(ev_ag)

#appropriations
ev_energy <- ergm(EVNETN~edges + nodematch("energycomm"), 
               control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_energy)
mcmc.diagnostics(ev_energy)
ev_energy_gof<-gof(ev_energy)
plot(ev_energy_gof)

#appropriations
ev_edwf <- ergm(EVNETN~edges + nodematch("educ_wf"), 
               control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_edwf)
mcmc.diagnostics(ev_edwf)
ev_edwf_gof<-gof(ev_edwf)
plot(ev_edwf)

# core 5
#appropriations
ev_core5 <- ergm(EVNETN~edges + 
                  nodematch("appropriations") + 
                  nodematch("help") + 
                  nodematch("agriculture") +
                  nodematch("energycomm") +
                  nodematch("educ_wf"), 
                control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_core5)
mcmc.diagnostics(ev_core5)
ev_core5_gof<-gof(ev_core5)
plot(ev_core5)



# STRUCTURAL                                                                                                              STRUCTURAL

# checking the degree distribution 
degreedist(EVNETN)
hist(degreedist(EVNETN))
# bimodal because of one-time policy experts and the those that sit in on panels very often 

# other distributions 
summary(EVNETN ~ edges + esp(1:40))
summary(EVNETN ~ edges + dsp(1:40))
summary(EVNETN ~ edges + nsp(1:40))

# ESP                                                                                                                      ESP

ev_esp10 <- ergm(EVNETN~edges + esp(0.10), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_esp10)
mcmc.diagnostics(ev_esp10)
ev_esp10_gof<-gof(ev_esp10)
plot(ev_esp10)

ev_esp25 <- ergm(EVNETN~edges + esp(0.25), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_esp25)
mcmc.diagnostics(ev_esp25)
ev_esp25_gof<-gof(ev_esp25)
plot(ev_esp25)

ev_esp50 <- ergm(EVNETN~edges + esp(0.50), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_esp50)
mcmc.diagnostics(ev_esp50)
ev_esp50_gof<-gof(ev_esp50)
plot(ev_esp50)

ev_esp75 <- ergm(EVNETN~edges + esp(0.75), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_esp75)
mcmc.diagnostics(ev_esp75)
ev_esp75_gof<-gof(ev_esp75)
plot(ev_esp75)


# GWESP                                                                                                                      GWESP

ev_gwesp10 <- ergm(EVNETN~edges + gwesp(0.10), 
                control=control.ergm(rivanna_cores, 
                MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwesp10)
mcmc.diagnostics(ev_gwesp10)
ev_gwesp10_gof<-gof(ev_gwesp10)
plot(ev_gwesp10)

ev_gwesp25 <- ergm(EVNETN~edges + gwesp(0.25), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwesp25)
mcmc.diagnostics(ev_gwesp25)
ev_gwesp25_gof<-gof(ev_gwesp25)
plot(ev_gwesp25)

ev_gwesp50 <- ergm(EVNETN~edges + gwesp(0.50), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwesp50)
mcmc.diagnostics(ev_gwesp50)
ev_gwesp50_gof<-gof(ev_gwesp50)
plot(ev_gwesp50)

ev_gwesp75 <- ergm(EVNETN~edges + gwesp(0.75), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwesp75)
mcmc.diagnostics(ev_gwesp75)
ev_gwesp75_gof<-gof(ev_gwesp75)
plot(ev_gwesp75)

# DSP                                                                                                                      DSP

ev_dsp10 <- ergm(EVNETN~edges + dsp(0.10), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_dsp10)
mcmc.diagnostics(ev_dsp10)
ev_dsp10_gof<-gof(ev_dsp10)
plot(ev_dsp10)

ev_dsp25 <- ergm(EVNETN~edges + dsp(0.25), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_dsp25)
mcmc.diagnostics(ev_dsp25)
ev_dsp25_gof<-gof(ev_dsp25)
plot(ev_dsp25)

ev_dsp50 <- ergm(EVNETN~edges + dsp(0.50), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_dsp50)
mcmc.diagnostics(ev_dsp50)
ev_dsp50_gof<-gof(ev_dsp50)
plot(ev_dsp50)

ev_dsp75 <- ergm(EVNETN~edges + dsp(0.75), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_dsp75)
mcmc.diagnostics(ev_dsp75)
ev_dsp75_gof<-gof(ev_dsp75)
plot(ev_dsp75)

# GWDSP                                                                                                                   GWDSP

ev_gwdsp10 <- ergm(EVNETN~edges + gwdsp(0.10), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwdsp10)
mcmc.diagnostics(ev_gwdsp10)
ev_gwdsp10_gof<-gof(ev_gwdsp10)
plot(ev_gwdsp10)

ev_gwdsp25 <- ergm(EVNETN~edges + gwdsp(0.25), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwdsp25)
mcmc.diagnostics(ev_gwdsp25)
ev_gwdsp25_gof<-gof(ev_gwdsp25)
plot(ev_gwdsp25)

ev_gwdsp50 <- ergm(EVNETN~edges + gwdsp(0.50), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwdsp50)
mcmc.diagnostics(ev_gwdsp50)
ev_gwdsp50_gof<-gof(ev_gwdsp50)
plot(ev_gwdsp50)

ev_gwdsp75 <- ergm(EVNETN~edges + gwdsp(0.75), 
                 control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev_gwdsp75)
mcmc.diagnostics(ev_gwdsp75)
ev_gwdsp75_gof<-gof(ev_gwdsp75)
plot(ev_gwdsp75)



















# https://github.com/statnet/Workshops/wiki
# https://statnet.org/Workshops/valued.html


#########################################################################
# Section 2: Final models (included in the manuscript)   
#########################################################################

# Final model 1: Looks at Structural Sweetspot of 4:11 degrees of connectivity 

ev.4_11d <- ergm(EVNETN~edges + degree(4:11), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d)

# Final model 2: Looks at best fitting model of nodal attributes 

ev.polfac.branchfac.core5 <- ergm(EVNETN~edges + nodefactor("political") + nodefactor("branch") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.polfac.branchfac.core5)

# Final model 3: Looks at structural sweetspot in combination with political party and House/Senate membership  

ev.4_11d.polfac.branchfac <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodefactor("branch"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.branchfac)

# Final model 4: Looks at structural sweetspot in combination with political party, branch membership, and impact of top-5 committees 

ev.4_11d.polfac.branchfac.core5 <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodefactor("branch") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.branchfac.core5)

# Final model summaries and GOF statistics 

final_m1 <- ev.4_11d
final_m2 <- ev.polfac.branchfac.core5
final_m3 <- ev.4_11d.polfac.branchfac
final_m4 <- ev.4_11d.polfac.branchfac.core5

summary(final_m1)
summary(final_m2)
summary(final_m3)
summary(final_m4)

gof_m1 <- gof(final_m1) 
gof_m2 <- gof(final_m2) 
gof_m3 <- gof(final_m3) 
gof_m4 <- gof(final_m4) 

plot(gof_m1); gof_m1
plot(gof_m2); gof_m2
plot(gof_m3); gof_m3
plot(gof_m4); gof_m4

# Interpreting the final models in terms of probability 

# Each number corresponds to the parameter included (order corresponds to model output summary). For example, including each parameters tells you how much the probability of a tie formating increases when accounting for 4-11 degrees of connectivity. Be mindful that it does not always make sense to include each parameter since actors are typically not assoicated with both the Democratic and Republication parties or in both the House and Senate (although there are some exceptions of both).

# This snippet tells your the combined coeffecient for the parameters in Model 1 (as well as which parameter aligns with which number)
exp(coef(final_m1)[c(1:9)])

# Model 1 with entire degree distribution 
plogis(sum(coef(final_m1)[c(1:9)])) 

# This looks at the probability of tie formation for a Democratic Senator on three of the five top committees
plogis(sum(coef(final_m2)[c(3, 4, 5, 6, 7)])) 

# Model 4 when policymakers have 4-6 degrees of connectivity and are: 

#(a) Democrats in the Senate and on two subcommittes 

plogis(sum(coef(final_m4)[c(1, 2, 3, 4, 10, 13, 15, 16)])) 

# (b) Republicans in the House and on two subcommittes 

plogis(sum(coef(final_m4)[c(1, 2, 3, 4, 11, 12, 17, 18)])) 

#########################################################################
# Section 3: Model development   
#########################################################################

# BASELINE MODEL 

ev <- ergm(EVNETN~edges, control=control.ergm(parallel=8, parallel.type="PSOCK"))
summary(ev)

ev$glm$aic
ev$glm
ev$theta1$coef
ev$mle.lik

# TESTING OF STRUCTURAL PARAMETERS 

# We started out by testing some structural features, including degree, weighted degree, k-star and alt-kstar. While I have included snippets of code for the first few models of each variable, more extensive testing was conducted on each of these variables as well as combinations of these variables before a more directed approach to variable selection was chosen based on degree distributions (see below).

# Testing degree ranges
ev.2deg <- ergm(EVNETN~edges + degree(2), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.2deg)

ev.3deg <- ergm(EVNETN~edges + degree(3), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.3deg)

ev.4deg <- ergm(EVNETN~edges + degree(4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.4deg)

ev.2_4deg <- ergm(EVNETN~edges + degree(2:4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.2_4deg)

# Testing weighted degree ranges
ev.2gwdeg <- ergm(EVNETN~edges + gwdegree(2), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.2gwdeg)

ev.3gwdeg <- ergm(EVNETN~edges + gwdegree(3), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.3gwdeg)

ev.4gwdeg <- ergm(EVNETN~edges + gwdegree(4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.4gwdeg)

ev.2_4gwdeg <- ergm(EVNETN~edges + gwdegree(2:4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.2_4gwdeg)

# Testing kstar parameters 

ev.k2 <- ergm(EVNETN~edges + kstar(2), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.k2)

ev.k3 <- ergm(EVNETN~edges + kstar(3), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.k3)

ev.k4 <- ergm(EVNETN~edges + kstar(4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.k4)

ev.k2_4 <- ergm(EVNETN~edges + kstar(2:4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.k4)

# Testing Altkstar

ev.altk2 <- ergm(EVNETN~edges + altkstar(2), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.k2_4)

ev.altk3 <- ergm(EVNETN~edges + altkstar(3), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.altk3)

ev.altk4 <- ergm(EVNETN~edges + altkstar(4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.altk4)

ev.altk2_4 <- ergm(EVNETN~edges + altkstar(2:4), control=control.ergm(rivanna_cores, parallel.type="PSOCK"))
summary(ev.altk2_4)

# Brandon ran several iterations of these models using a variety of ranges and combinations of structural parameters, including degree, gwdegree (both measures of structural connectivity), kstar, altkstar (both proxies for brokering), isolates, twopath, meandeg, balance, and sociality. While most of these mechanisms were of theoretical interest (especially the degree and brokering parameters), few of these parameters ultimately generated models that converged except for models predicated on degree distributions. The biggest issue here is that we assumed based on the general network summary statistics that the sparseness of the network would mean low degree ranges would be our best predictors for the model. Eventually, we were able to get a better handle on which degree ranges were actually most appropriate by looking at the summary statistics of the degree ranges. 

summary(EVNETN ~ edges + degree(1:21))

# The output demonstrates that degree ranges 3 through 14 are the most likely to provide the best fitting model based on degree distributions. However, as you will see more below, even after we established this range of connectivity, we still have to test various iterations of the degree parameters to find the best fitting models. In the end, we found that the range was actually 4:11 actors - point to a "structural sweetspot" where actors must be at least connected to four other actors to increase their probability of sharing evidence, but do not necessarily need to be connected to more than 11 to raise that probability. 

# TESTING ATTRIBUTE-BASED MODELS  

# Next, we began testing various how various nodal attributes affected the probability of tie formation. Intitially, these models were predicated on affiliation (affil), a proxy for whether an actor (polnotpol), 
#  FINISH THIS HERE !!!!!!! 
# Note that some of these models were developed before the degree distribution was honed more precisely. However, at the time of development, it was clear that four degrees of connection would have a substantial impact on the model. 

# Affiliation models 

ev.affil <- ergm(EVNETN~edges + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.affil)
mcmc.diagnostics(ev.affil)
evgof.1.affil<-gof(ev.affil)
plot(ev.affil)

ev.d4.bal.affil <- ergm(EVNETN~edges + degree(2) + balance + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.bal.affil)
mcmc.diagnostics(ev.d4.bal.affil)
evgof.1.d4.bal.affil<-gof(ev.d4.bal.affil)
plot(evgof.1.d4.bal.affil)

ev.d4.affil <- ergm(EVNETN~edges + degree(4) + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.affil)
mcmc.diagnostics(ev.d4.affil)
evgof.1.d4.affil<-gof(ev.d4.affil)
plot(evgof.1.d4.affil)

ev.k4.affil <- ergm(EVNETN~edges + kstar(4) + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.k4.affil)
mcmc.diagnostics(ev.k4.affil)
evgof.1.k4.affil<-gof(ev.k4.affil)
plot(evgof.1.k4.affil)

# Gender models 

ev.d4.gend <- ergm(EVNETN~edges + degree(4) + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.gend)
mcmc.diagnostics(ev.d4.gend)
evgof.1.d4.gend<-gof(ev.d4.gend)
plot(evgof.1.d4.gend)

ev.k4.gend <- ergm(EVNETN~edges + kstar(4) + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.k4.gend)
mcmc.diagnostics(ev.k4.gend)
evgof.1.k4.gend<-gof(ev.k4.gend)
plot(evgof.1.k4.gend)

ev.iso.gen <- ergm(EVNETN~edges + isolates + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.iso.gen)
mcmc.diagnostics(ev.iso.gen)
evgof.1.iso.gen<-gof(ev.iso.gen)
plot(evgof.1.iso.gen)

# Tenure models 

ev.d4.tenure <- ergm(EVNETN~edges + degree(4) + nodematch("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.tenure)
mcmc.diagnostics(ev.d4.tenure)
evgof.1.d4.tenure<-gof(ev.d4.tenure)
plot(evgof.1.d4.tenure)

ev.k4.tenure <- ergm(EVNETN~edges + kstar(4) + nodematch("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.k4.tenure)
mcmc.diagnostics(ev.k4.tenure)
evgof.1.k4.tenure<-gof(ev.k4.tenure)
plot(evgof.1.k4.tenure)

ev.d4.k4.tenure <- ergm(EVNETN~edges + degree(4) + kstar(4) + nodematch("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.k4.tenure)
mcmc.diagnostics(ev.d4.k4.tenure)
evgof.1.d4.k4.tenure<-gof(ev.d4.k4.tenure)
plot(evgof.1.d4.k4.tenure)

# Political proxy models 

ev.pnp <- ergm(EVNETN~edges + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.pnp)
mcmc.diagnostics(ev.pnp)
evgof.1.pnp<-gof(ev.pnp)
plot(evgof.1.pnp)

ev.d4.pnp <- ergm(EVNETN~edges + degree(4) + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.pnp)
mcmc.diagnostics(ev.d4.pnp)
evgof.1.d4.pnp<-gof(ev.d4.pnp)
plot(evgof.1.d4.pnp)

ev.k4.pnp <- ergm(EVNETN~edges + kstar(4) + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.k4.pnp)
mcmc.diagnostics(ev.k4.pnp)
evgof.1.k4.pnp<-gof(ev.k4.pnp)
plot(evgof.1.k4.pnp)

ev.k4.iso.pnp <- ergm(EVNETN~edges + kstar(4) + isolates + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.k4.iso.pnp)
mcmc.diagnostics(ev.k4.iso.pnp)
evgof.1.k4.iso.pnp<-gof(ev.k4.iso.pnp)
plot(evgof.1.k4.iso.pnp)

ev.d4.k4.pnp <- ergm(EVNETN~edges + degree(4) + kstar(4) + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.k4.pnp)
mcmc.diagnostics(ev.d4.k4.pnp)
evgof.1.d4.k4.pnp<-gof(ev.d4.k4.pnp)
plot(evgof.1.d4.k4.pnp)

# Permutations of nodal attributes with four degrees of connectivity 

ev.d4.affil.gend <- ergm(EVNETN~edges + degree(4) + nodematch("affil") + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.affil.gend)
mcmc.diagnostics(ev.d4.affil.gend)
evgof.1.d4.affil.gend<-gof(ev.d4.affil.gend)
plot(evgof.1.d4.affil.gend)

ev.d4.affil.pol <- ergm(EVNETN~edges + degree(4) + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.affil.pol)
mcmc.diagnostics(ev.d4.affil.pol)
evgof.1.d4.affil.pol<-gof(ev.d4.affil.pol)
plot(evgof.1.d4.affil.pol)

ev.d4.affil.gend.pol <- ergm(EVNETN~edges + degree(4) + nodematch("affil") + nodematch("gender") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.affil.gend.pol)
mcmc.diagnostics(ev.d4.affil.gend.pol)
evgof.1.d4.affil.gend.pol<-gof(ev.d4.affil.gend.pol)
plot(evgof.1.d4.affil.gend.pol)

ev.d4.affil.gend.ten <- ergm(EVNETN~edges + degree(4) + nodematch("affil") + nodematch("gender") + nodematch("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.affil.gend.ten)
mcmc.diagnostics(ev.d4.affil.gend.ten)
evgof.1.d4.affil.gend.ten<-gof(ev.d4.affil.gend.ten)
plot(evgof.1.d4.affil.gend.ten)

ev.d4.affil.gend.pol.ten <- ergm(EVNETN~edges + degree(4) + nodematch("affil") + nodematch("gender") + nodematch("political") + nodematch("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.affil.gend.pol.ten)
mcmc.diagnostics(ev.d4.affil.gend.pol.ten)
evgof.1.d4.affil.gend.pol.ten<-gof(ev.d4.affil.gend.pol.ten)
plot(evgof.1.d4.affil.gend.pol.ten)

ev.d4.covten <- ergm(EVNETN~edges + degree(4) + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.covten)
mcmc.diagnostics(ev.d4.covten)
evgof.1.d4.covten<-gof(ev.d4.covten)
plot(evgof.1.d4.covten)

ev.k4.covten <- ergm(EVNETN~edges + kstar(4) + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.k4.covten)
mcmc.diagnostics(ev.k4.covten)
evgof.1.k4.covten<-gof(ev.k4.covten)
plot(evgof.1.k4.covten)

ev.d4.k4.covten <- ergm(EVNETN~edges + degree(4) + kstar(4) + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.d4.k4.covten)
mcmc.diagnostics(ev.d4.k4.covten)
evgof.1.d4.k4.covten<-gof(ev.d4.k4.covten)
plot(evgof.1.d4.k4.covten)

ev.4d.covten.pnp <- ergm(EVNETN~edges + degree(4) + nodecov("tenure") + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4d.covten.pnp)
mcmc.diagnostics(ev.4d.covten.pnp)
evgof.1.4d.covten.pnp<-gof(ev.4d.covten.pnp)
plot(evgof.1.4d.covten.pnp)

ev.4d.md.covten.pnp <- ergm(EVNETN~edges + degree(4) + meandeg + nodecov("tenure") + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4d.md.covten.pnp)
mcmc.diagnostics(ev.4d.md.covten.pnp)
evgof.1.4d.md.covten.pnp<-gof(ev.4d.md.covten.pnp)
plot(evgof.1.4d.md.covten.pnp)

# STRUCTURAL + NODAL ATTRIBUTE-BASED MODELS 

# As mentioned above, Brandon's development of the model assumed that low degree parameters would be the best fit because the network was fairly sparse. However, after running more descriptives on the network and honing the nodal attributes further, we realized that a degree distribution would be a better fit. Below, I ran various iterations of these degree parameters. 

deg_dist <- summary(EVNETN ~ edges + degree(1:25))
deg_dist <- c(6,3,20,25,12,19,4,20,30,13,35,14, 9,10,4,4,5,3,3,17,2 ,5,1,2,4)
hist(deg_dist)

ev.3_6d.20d.covten.aff.pol.gend <- ergm(EVNETN~edges + degree(4:6) + degree(20) + nodecov("tenure")  + nodematch("affil") + nodematch("political") + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_6d.20d.covten.aff.pol.gend)
mcmc.diagnostics(ev.3_6d.20d.covten.aff.pol.gend)
evgof.1.3_6d.20d.covten.aff.pol.gend<-gof(ev.3_6d.20d.covten.aff.pol.gend)
plot(evgof.1.3_6d.20d.covten.aff.pol.gend)

ev.3_6d.8_12d.20d.covten.aff.pol.gend <- ergm(EVNETN~edges + degree(3:6) + degree(8:12) + degree(20) + nodecov("tenure")  + nodematch("affil") + nodematch("political") + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_6d.8_12d.20d.covten.aff.pol.gend)
mcmc.diagnostics(ev.3_6d.8_12d.20d.covten.aff.pol.gend)
evgof.1.3_6d.8_12d.20d.covten.aff.pol.gend<-gof(ev.3_6d.8_12d.20d.covten.aff.pol.gend)
plot(evgof.1.3_6d.8_12d.20d.covten.aff.pol.gend)

ev.8_12d.covten.aff.pol.gend <- ergm(EVNETN~edges + degree(8:12) + nodecov("tenure")  + nodematch("affil") + nodematch("political") + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.8_12d.covten.aff.pol.gend)
mcmc.diagnostics(ev.4_6d.8_12d.20d.covten.aff.pol.gend)
evgof.1.4_6d.8_12d.20d.covten.aff.pol.gend<-gof(ev.4_6d.8_12d.20d.covten.aff.pol.gend)
plot(evgof.1.4_6d.8_12d.20d.covten.aff.pol.gend)

ev.5_6d.8_12d.20d.covten.aff.pol.gend <- ergm(EVNETN~edges + degree(5:6) + degree(8:12) + degree(20) + nodecov("tenure")  + nodematch("affil") + nodematch("political") + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.5_6d.8_12d.20d.covten.aff.pol.gend)
mcmc.diagnostics(ev.5_6d.8_12d.20d.covten.aff.pol.gend)
evgof.1.5_6d.8_12d.20d.covten.aff.pol.gend<-gof(ev.5_6d.8_12d.20d.covten.aff.pol.gend)
plot(evgof.1.5_6d.8_12d.20d.covten.aff.pol.gend)

ev.4_6d.8_11d.20d.covten.aff.pol.gend <- ergm(EVNETN~edges + degree(4:6) + degree(8:11) + degree(20) + nodecov("tenure")  + nodematch("affil") + nodematch("political") + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_6d.8_11d.20d.covten.aff.pol.gend)
mcmc.diagnostics(ev.4_6d.8_11d.20d.covten.aff.pol.gend)
evgof.1.4_6d.8_11d.20d.covten.aff.pol.gend<-gof(ev.4_6d.8_11d.20d.covten.aff.pol.gend)
plot(evgof.1.4_6d.8_11d.20d.covten.aff.pol.gend)

# What if I just include the terms with 20 or more degree ties?

ev.3_4d.6d.8_9d.11d.20d.covten.aff.pol.gend <- ergm(EVNETN~edges + degree(3:4) + degree(6) + degree(8:9) + degree(11) + degree(20) + nodecov("tenure")  + nodematch("affil") + nodematch("political") + nodematch("gender"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_4d.6d.8_9d.11d.20d.covten.aff.pol.gend)
mcmc.diagnostics(ev.3_4d.6d.8_9d.11d.20d.covten.aff.pol.gend)
evgof.1.3_4d.6d.8_9d.11d.20d.covten.aff.pol.gend<-gof(ev.3_4d.6d.8_9d.11d.20d.covten.aff.pol.gend)
plot(evgof.1.3_4d.6d.8_9d.11d.20d.covten.aff.pol.gend)

# Trying to fit degree(3:20)

ev.3_20d.covten.aff <- ergm(EVNETN~edges + degree(3:20) + nodecov("tenure")  + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_20d.covten.aff)
mcmc.diagnostics(ev.3_20d.covten.aff)
evgof.1.3_20d.covten.aff<-gof(ev.3_20d.covten.aff)
plot(evgof.1.3_20d.covten.aff)

# Trying to fit degree(3:12)

ev.3_12d.covten.aff <- ergm(EVNETN~edges + degree(3:12) + nodecov("tenure")  + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_12d.covten.aff)
mcmc.diagnostics(ev.3_12d.covten.aff)
evgof.1.3_12d.covten.aff<-gof(ev.3_12d.covten.aff)
plot(evgof.1.3_12d.covten.aff)

ev.3_12d.covten.aff.pol <- ergm(EVNETN~edges + degree(3:12) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_12d.covten.aff.pol)
mcmc.diagnostics(ev.3_12d.covten.aff.pol)
evgof.1.3_12d.covten.aff.pol<-gof(ev.3_12d.covten.aff.pol)
plot(evgof.1.3_12d.covten.aff.pol)

ev.3_11d.covten.aff.pol <- ergm(EVNETN~edges + degree(3:11) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.covten.aff.pol)
mcmc.diagnostics(ev.3_11d.covten.aff.pol)
evgof.1.3_11d.covten.aff.pol<-gof(ev.3_11d.covten.aff.pol)
plot(evgof.1.3_11d.covten.aff.pol)

ev.3_11d.covten.aff <- ergm(EVNETN~edges + degree(3:11) + nodecov("tenure")  + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.covten.aff)
mcmc.diagnostics(ev.3_11d.covten.aff)
evgof.1.3_11d.covten.aff<-gof(ev.3_11d.covten.aff)
plot(evgof.1.3_11d.covten.aff)

# Looks like 3:11 is going to be pretty close to the final distribution so I continued to hone the attribute parameters... 

# nodefactor terms 

ev.3_11d.nf_aff.ncov_ten.pol <- ergm(EVNETN~edges + degree(3:11) + nodefactor("affil") + nodecov("tenure") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.nf_aff.ncov_ten.pol)
mcmc.diagnostics(ev.3_11d.nf_aff.ncov_ten.pol)
evgof.1.3_11d.nf_aff.ncov_ten.pol<-gof(ev.3_11d.nf_aff.ncov_ten.pol)
plot(evgof.1.3_11d.nf_aff.ncov_ten.pol)

ev.3_11d.aff.ncov_ten.nf_pol <- ergm(EVNETN~edges + degree(3:11) + nodematch("affil") + nodecov("tenure") + nodefactor("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.aff.ncov_ten.nf_pol)
mcmc.diagnostics(ev.3_11d.aff.ncov_ten.nf_pol)
evgof.1.3_11d.aff.ncov_ten.nf_pol<-gof(ev.3_11d.aff.ncov_ten.nf_pol)
plot(evgof.1.3_11d.aff.ncov_ten.nf_pol)

ev.3_11d.nf_aff.ncov_ten.nf_pol <- ergm(EVNETN~edges + degree(3:11) + nodefactor("affil") + nodecov("tenure") + nodefactor("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.nf_aff.ncov_ten.nf_pol)
mcmc.diagnostics(ev.3_11d.nf_aff.ncov_ten.nf_pol)
evgof.1.3_11d.nf_aff.ncov_ten.nf_pol<-gof(ev.3_11d.nf_aff.ncov_ten.nf_pol)
plot(evgof.1.3_11d.nf_aff.ncov_ten.nf_pol)

# nodemix terms 

ev.3_11d.nm_affil.pol.nc_ten <- ergm(EVNETN~edges + degree(3:11) + nodemix("affil", base = 0) + nodematch("political") + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.nm_affil.pol.nc_ten)
mcmc.diagnostics(ev.3_11d.nm_affil.pol.nc_ten)
evgof.1.3_11d.nm_affil.pol.nc_ten<-gof(ev.3_11d.nm_affil.pol.nc_ten)
plot(evgof.1.3_11d.nm_affil.pol.nc_ten)

ev.3_11d.nm_pol.affil.nc_ten <- ergm(EVNETN~edges + degree(3:11) + nodemix("political", base = 0) + nodematch("affil") + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.nm_pol.affil.nc_ten)
mcmc.diagnostics(ev.3_11d.nm_pol.affil.nc_ten)
evgof.1.3_11d.nm_pol.affil.nc_ten<-gof(ev.3_11d.nm_pol.affil.nc_ten)
plot(evgof.1.3_11d.nm_pol.affil.nc_ten)

ev.3_11d.nm_pol.nm_affil.nc_ten <- ergm(EVNETN~edges + degree(3:11) + nodemix("affil", base = 0) + nodemix("political", base=0) + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.nm_pol.nm_affil.nc_ten)
mcmc.diagnostics(ev.3_11d.nm_pol.nm_affil.nc_ten)
evgof.1.3_11d.nm_pol.nm_affil.nc_ten<-gof(ev.3_11d.nm_pol.nm_affil.nc_ten)
plot(evgof.1.3_11d.nm_pol.nm_affil.nc_ten)

# absdiff terms 

ev.3_11d.ad_affil.pol.nc_ten <- ergm(EVNETN~edges + degree(3:11) + absdiff("affil") + nodematch("political") + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.ad_affil.pol.nc_ten)
mcmc.diagnostics(ev.3_11d.ad_affil.pol.nc_ten)
evgof.1.3_11d.ad_affil.pol.nc_ten<-gof(ev.3_11d.ad_affil.pol.nc_ten)
plot(evgof.1.3_11d.ad_affil.pol.nc_ten)

ev.3_11d.ad_pol.affil.nc_ten <- ergm(EVNETN~edges + degree(3:11) + absdiff("political") + nodematch("affil") + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.ad_pol.affil.nc_ten)
mcmc.diagnostics(ev.3_11d.ad_pol.affil.nc_ten)
evgof.1.3_11d.ad_pol.affil.nc_ten<-gof(ev.3_11d.ad_pol.affil.nc_ten)
plot(evgof.1.3_11d.ad_pol.affil.nc_ten)

ev.3_11d.ad_pol.ad_affil.nc_ten <- ergm(EVNETN~edges + degree(3:11) + absdiff("political") + absdiff("affil") + nodecov("tenure"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.ad_pol.ad_affil.nc_ten)
mcmc.diagnostics(ev.3_11d.ad_pol.ad_affil.nc_ten)
evgof.1.3_11d.ad_pol.ad_affil.nc_ten<-gof(ev.3_11d.ad_pol.ad_affil.nc_ten)
plot(evgof.1.3_11d.ad_pol.ad_affil.nc_ten)

# At this point, it looked like 3:11 would be the best fit in terms of the degree distribution, which is only slightly off of the final parameters used in the manuscript (i.e. 4:11). Given the success of this approach, Brandon wanted to re-try some of the models with a descriptives based approach to improving the GOF. After running descriptives for edgewise-shared partners (both weighted and non-weighted), dyad-shared partners and non-edgewise shared partners, we developed several models that tested the esp, gwesp, nsp and dsp terms. 

#summary of degree tests: either  the 4:11 or 6:11 models should fit best 

summary(EVNETN ~ edges + degree(1:25))
summary(EVNETN ~ edges + esp(1:20))
summary(EVNETN ~ edges + dsp(1:20))
summary(EVNETN ~ edges + nsp(1:20))
summary(EVNETN ~ edges + graphletCount(1:29)) #6 & 18 = brokerage in Lusher, #9:11 = brokerage for Yaveroglu (both are correct)

ev.gwesp25.ten.aff.pol <- ergm(EVNETN~edges + gwesp(0.25,fixed=T) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.gwesp25.ten.aff.pol)
mcmc.diagnostics(ev.gwesp25.ten.aff.pol)
evgof.1.gwesp25.ten.aff.pol<-gof(ev.gwesp25.ten.aff.pol)
plot(evgof.1.gwesp25.ten.aff.pol)

ev.g50.3_11d.covten.aff.pol <- ergm(EVNETN~edges + gwesp(0.50,fixed=T) + degree(3:11) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.g50.3_11d.covten.aff.pol)
mcmc.diagnostics(ev.g50.3_11d.covten.aff.pol)
evgof.1.g50.3_11d.covten.aff.pol<-gof(ev.g50.3_11d.covten.aff.pol)
plot(evgof.1.g50.3_11d.covten.aff.pol)

ev.g75.3_11d.covten.aff.pol <- ergm(EVNETN~edges + gwesp(0.75,fixed=T) + degree(3:11) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.g75.3_11d.covten.aff.pol)
mcmc.diagnostics(ev.g75.3_11d.covten.aff.pol)
evgof.1.g75.3_11d.covten.aff.pol<-gof(ev.g75.3_11d.covten.aff.pol)
plot(evgof.1.g75.3_11d.covten.aff.pol)

ev.esp10.ten.aff.pol <- ergm(EVNETN~edges + esp(0.10) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.esp10.ten.aff.pol)
mcmc.diagnostics(ev.esp10.ten.aff.pol)
evgof.1.esp10.ten.aff.pol<-gof(ev.esp10.ten.aff.pol)
plot(evgof.1.esp10.ten.aff.pol)

ev.esp25.ten.aff.pol <- ergm(EVNETN~edges + esp(0.25) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.esp25.ten.aff.pol)
mcmc.diagnostics(ev.esp25.ten.aff.pol)
evgof.1.esp25.ten.aff.pol<-gof(ev.esp25.ten.aff.pol)
plot(evgof.1.esp25.ten.aff.pol)

ev.esp75.3_11d.ten.aff.pol <- ergm(EVNETN~edges + esp(0.75) + degree(3:11) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.esp75.3_11d.ten.aff.pol)
mcmc.diagnostics(ev.esp75.3_11d.ten.aff.pol)
evgof.1.esp75.3_11d.ten.aff.pol<-gof(ev.esp75.3_11d.ten.aff.pol)
plot(evgof.1.esp75.3_11d.ten.aff.pol)

ev.esp25.gwdsp25.nc_ten.aff.pol <- ergm(EVNETN~edges + esp(0.25) + gwdsp(0.25, fixed = FALSE) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.esp25.gwdsp25.nc_ten.aff.pol)
mcmc.diagnostics(ev.esp25.gwdsp25.nc_ten.aff.pol)
evgof.1.esp25.gwdsp25.nc_ten.aff.pol<-gof(ev.esp25.gwdsp25.nc_ten.aff.pol)
plot(evgof.1.esp25.gwdsp25.nc_ten.aff.pol)

ev.gwdsp25.gwesp25.nc_ten.aff.pol <- ergm(EVNETN~edges + gwdsp(0.25, fixed = FALSE) + gwesp(0.25, fixed = FALSE) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.gwdsp25.gwesp25.nc_ten.aff.pol)
mcmc.diagnostics(ev.gwdsp25.gwesp25.nc_ten.aff.pol)
evgof.1.gwdsp25.gwesp25.nc_ten.aff.pol<-gof(ev.gwdsp25.gwesp25.nc_ten.aff.pol)
plot(evgof.1.gwdsp25.gwesp25.nc_ten.aff.pol)

ev.gwdsp25.3_11d.nc_ten.aff.pol <- ergm(EVNETN~edges + gwdsp(0.25, fixed = FALSE) + degree(3:11) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.gwdsp25.3_11d.nc_ten.aff.pol)
mcmc.diagnostics(ev.gwdsp25.3_11d.nc_ten.aff.pol)
evgof.1.gwdsp25.3_11d.nc_ten.aff.pol<-gof(ev.gwdsp25.3_11d.nc_ten.aff.pol)
plot(evgof.1.gwdsp25.3_11d.nc_ten.aff.pol)

ev.3_11d.covten.aff.pol <- ergm(EVNETN~edges + degree(3:11) + nodecov("tenure") + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.3_11d.covten.aff.pol)
mcmc.diagnostics(ev.3_11d.covten.aff.pol)
evgof.3_11d.covten.aff.pol<-gof(ev.3_11d.covten.aff.pol)
plot(evgof.3_11d.covten.aff.pol)
evgof.3_11d.covten.aff.pol

ev.esp25.covten.aff.pol <- ergm(EVNETN~edges + esp(0.25) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.esp25.covten.aff.pol)
mcmc.diagnostics(ev.esp25.covten.aff.pol)
evgof.esp25.covten.aff.pol<-gof(ev.esp25.covten.aff.pol)
plot(evgof.esp25.covten.aff.pol)
evgof.esp25.covten.aff.pol

ev.gwesp50.3_11d.covten.aff.pol <- ergm(EVNETN~edges + gwesp(0.50) + degree(3:11) + nodecov("tenure")  + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.gwesp50.3_11d.covten.aff.pol)
mcmc.diagnostics(ev.gwesp50.3_11d.covten.aff.pol)
evgof.gwesp50.3_11d.covten.aff.pol<-gof(ev.gwesp50.3_11d.covten.aff.pol)
plot(evgof.gwesp50.3_11d.covten.aff.pol)
evgof.gwesp50.3_11d.covten.aff.pol

ev.esp25.dsp3 <- ergm(EVNETN~edges + esp(0.25) + dsp(3), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.esp25.dsp3)
mcmc.diagnostics(ev.esp25.dsp3)
evgof.esp25.dsp3<-gof(ev.esp25.dsp3)
plot(evgof.esp25.dsp3)
evgof.esp25.dsp3

ev.esp25.dsp10 <- ergm(EVNETN~edges + esp(0.25) + dsp(10), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.esp25.dsp10)
mcmc.diagnostics(ev.esp25.dsp10)
evgof.esp25.dsp10<-gof(ev.esp25.dsp10)
plot(evgof.esp25.dsp10)
evgof.esp25.dsp10

# While the gwesp, nsp and dsp models did not work out well (most didn't even converge), the lower esp-based models did proffer some interesting results so I went back and developed these models a bit further with the degree distributions inserted into the models. 

summary(EVNETN ~ edges + degree(1:21))

ev.8_11d.covten.aff.pol <- ergm(EVNETN~edges + degree(8:11) + nodecov("tenure") + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.8_11d.covten.aff.pol)
evgof.8_11d.covten.aff.pol<-gof(ev.8_11d.covten.aff.pol)
plot(evgof.8_11d.covten.aff.pol); evgof.8_11d.covten.aff.pol

ev.6_11d.covten.aff.pol <- ergm(EVNETN~edges + degree(6:11) + nodecov("tenure") + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.6_11d.covten.aff.pol)
evgof.6_11d.covten.aff.pol<-gof(ev.6_11d.covten.aff.pol)
plot(evgof.6_11d.covten.aff.pol); evgof.6_11d.covten.aff.pol

ev.4_11d.covten.aff.pol <- ergm(EVNETN~edges + degree(4:11) + nodecov("tenure") + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.covten.aff.pol)
evgof.4_11d.covten.aff.pol<-gof(ev.4_11d.covten.aff.pol)
plot(evgof.4_11d.covten.aff.pol); evgof.4_11d.covten.aff.pol

ev.4_11d.esp10.covten.aff.pol <- ergm(EVNETN~edges + degree(4:11) + esp(10) + nodecov("tenure") + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.esp10.covten.aff.pol)
evgof.4_11d.esp10.covten.aff.pol<-gof(ev.4_11d.esp10.covten.aff.pol)
plot(evgof.4_11d.esp10.covten.aff.pol); evgof.4_11d.esp10.covten.aff.pol

ev.4_1120d.esp7_1320.covten.aff.pol <- ergm(EVNETN~edges + degree(4:11) + degree(20) + esp(7:13) + esp(19) + nodecov("tenure") + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_1120d.esp7_1320.covten.aff.pol)
evgof.4_1120d.esp7_1320.covten.aff.pol<-gof(ev.4_1120d.esp7_1320.covten.aff.pol)
plot(evgof.4_1120d.esp7_1320.covten.aff.pol); evgof.4_1120d.esp7_1320.covten.aff.pol

ev.4_1120d.esp3_1320.nsp1_2.aff.pol <- ergm(EVNETN~edges + degree(4:11) + degree(20) + esp(3:13) + esp(19) + nsp(1:2) + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_1120d.esp3_1320.nsp1_2.aff.pol)
evgof.4_1120d.esp3_1320.nsp1_2.aff.pol<-gof(ev.4_1120d.esp3_1320.nsp1_2.aff.pol)
plot(evgof.4_1120d.esp3_1320.nsp1_2.aff.pol); evgof.4_1120d.esp3_1320.nsp1_2.aff.pol

ev.4_1120d.esp5_1320.nsp1_2.dsp1_2.aff.pol <- ergm(EVNETN~edges + degree(4:11) + degree(20) + esp(5:13) + esp(19) + nsp(1:2) + dsp(1:2) + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_1120d.esp5_1320.nsp1_2.dsp1_2.aff.pol)
evgof.4_1120d.esp5_1320.nsp1_2.dsp1_2.aff.pol<-gof(ev.4_1120d.esp5_1320.nsp1_2.dsp1_2.aff.pol)
plot(evgof.4_1120d.esp5_1320.nsp1_2.dsp1_2.aff.pol); evgof.4_1120d.esp5_1320.nsp1_2.dsp1_2.aff.pol

ev.14_111420d.esp5_1320.nsp1_2.dsp1_2 <- ergm(EVNETN~edges + degree(1) + degree(4:11) + degree(14) + degree(20) + esp(5:13) + esp(19) + nsp(1:2) + dsp(1:2), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.14_111420d.esp5_1320.nsp1_2.dsp1_2)
evgof.14_111420d.esp5_1320.nsp1_2.dsp1_2<-gof(ev.14_111420d.esp5_1320.nsp1_2.dsp1_2)
plot(evgof.14_111420d.esp5_1320.nsp1_2.dsp1_2); evgof.14_111420d.esp5_1320.nsp1_2.dsp1_2

ev.4_1120d.esp3_1320.nsp1_2.dsp1_2.aff.pol <- ergm(EVNETN~edges + degree(4:11) + degree(20) + esp(3:13) + esp(19) + nsp(1:2) + dsp(1:2) + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_1120d.esp3_1320.nsp1_2.dsp1_2.aff.pol)
evgof.4_1120d.esp3_1320.nsp1_2.dsp1_2.aff.pol<-gof(ev.4_1120d.esp3_1320.nsp1_2.dsp1_2.aff.pol)
plot(evgof.4_1120d.esp3_1320.nsp1_2.dsp1_2.aff.pol); evgof.4_1120d.esp3_1320.nsp1_2.dsp1_2.aff.pol

ev.4_11d.nsp1_2.dsp1_2.aff.pol <- ergm(EVNETN~edges + degree(4:11) + nsp(1:2) + dsp(1:2) + nodematch("affil") + nodematch("political"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.nsp1_2.dsp1_2.aff.pol)
evgof.4_11d.nsp1_2.dsp1_2.aff.pol<-gof(ev.4_11d.nsp1_2.dsp1_2.aff.pol)
plot(evgof.4_11d.nsp1_2.dsp1_2.aff.pol); evgof.4_11d.nsp1_2.dsp1_2.aff.pol

# Eventually, I decided just to remove all of the nodal attributes and try to develop a model solely based on structural features. 

ev.14_111420d.esp5_1320.nsp1_2.dsp1_2 <- ergm(EVNETN~edges + degree(1) + degree(4:11) + degree(14) + degree(20) + esp(5:13) + esp(19) + nsp(1:2) + dsp(1:2), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.14_111420d.esp5_1320.nsp1_2.dsp1_2)
evgof.14_111420d.esp5_1320.nsp1_2.dsp1_2<-gof(ev.14_111420d.esp5_1320.nsp1_2.dsp1_2)
plot(evgof.14_111420d.esp5_1320.nsp1_2.dsp1_2); evgof.14_111420d.esp5_1320.nsp1_2.dsp1_2

ev.4_1120d.esp5_1320.nsp1.dsp1 <- ergm(EVNETN~edges + degree(4:11) + degree(20) + esp(5:13) + esp(19) + nsp(1) + dsp(1), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_1120d.esp5_1320.nsp1.dsp1)
evgof.4_1120d.esp5_1320.nsp1.dsp1<-gof(ev.4_1120d.esp5_1320.nsp1.dsp1)
plot(evgof.4_1120d.esp5_1320.nsp1.dsp1); evgof.4_1120d.esp5_1320.nsp1.dsp1

ev.4_6d.8_12d.20d.esp5_1319.pol.pnp <- ergm(EVNETN~edges + degree(4:6) + degree(8:12) + degree(20) + esp(5:13) + esp(19) + nodematch("political") + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_6d.8_12d.20d.esp5_1319.pol.pnp)
evgof.4_6d.8_12d.20d.esp5_1319.pol.pnp<-gof(ev.4_6d.8_12d.20d.esp5_1319.pol.pnp)
plot(evgof.4_6d.8_12d.20d.esp5_1319.pol.pnp); evgof.4_6d.8_12d.20d.esp5_1319.pol.pnp

ev.4_6d.8_12d.20d.esp5_1319.pol.pnp.core5 <- ergm(EVNETN~edges + degree(4:6) + degree(8:12) + degree(20) + esp(5:13) + esp(19) + nodematch("political") + nodematch("polnotpol") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_6d.8_12d.20d.esp5_1319.pol.pnp.core5)
evgof.4_6d.8_12d.20d.esp5_1319.pol.pnp.core5<-gof(ev.4_6d.8_12d.20d.esp5_1319.pol.pnp.core5)
plot(evgof.4_6d.8_12d.20d.esp5_1319.pol.pnp.core5); evgof.4_6d.8_12d.20d.esp5_1319.pol.pnp.core5

ev.4_6d.8_12d.20d.esp5_1319.aff <- ergm(EVNETN~edges + degree(4:6) + degree(8:12) + degree(20) + esp(5:13) + esp(19) + nodematch("affil"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_6d.8_12d.20d.esp5_1319.aff)
evgof.4_6d.8_12d.20d.esp5_1319.aff<-gof(ev.4_6d.8_12d.20d.esp5_1319.aff)
plot(evgof.4_6d.8_12d.20d.esp5_1319.aff); evgof.4_6d.8_12d.20d.esp5_1319.aff

ev.4_6d.8_12d.20d.esp5_1319.aff.core5 <- ergm(EVNETN~edges + degree(4:6) + degree(8:12) + degree(20) + esp(5:13) + esp(19) + nodematch("affil") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_6d.8_12d.20d.esp5_1319.aff.core5)
evgof.4_6d.8_12d.20d.esp5_1319.aff.core5<-gof(ev.4_6d.8_12d.20d.esp5_1319.aff.core5)
plot(evgof.4_6d.8_12d.20d.esp5_1319.aff.core5); evgof.4_6d.8_12d.20d.esp5_1319.aff.core5

# This model is solely based on structural parameters and provides a really low AIC/BIC with . The biggest problem is that this model (1) produces an overfit model ---- and (b) that it does not provide any clear policy implications. What can we make of a model that says structural factors are the only parameters that influence evidence sharing in policymaking contexts? 

ev.4_6d.8_12d.20d.esp5_1319 <- ergm(EVNETN~edges + degree(4:6) + degree(8:12) + degree(20) + esp(5:13) + esp(19), control=control.ergm(parallel=12, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_6d.8_12d.20d.esp5_1319)
evgof.4_6d.8_12d.20d.esp5_1319<-gof(ev.4_6d.8_12d.20d.esp5_1319)
plot(evgof.4_6d.8_12d.20d.esp5_1319); evgof.4_6d.8_12d.20d.esp5_1319

ev.4_11d.esp5_13 <- ergm(EVNETN~edges + degree(4:11) + esp(5:13), control=control.ergm(parallel=12, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.esp5_13)
evgof.4_11d.esp5_13<-gof(ev.4_11d.esp5_13)
plot(evgof.4_11d.esp5_13); evgof.4_11d.esp5_13

ev.4_11d.esp5_13.polfac.branchfac <- ergm(EVNETN~edges + degree(4:11) + esp(5:13) + nodefactor("political") + nodefactor("branch"), control=control.ergm(parallel=12, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.esp5_13.polfac.branchfac)
evgof.4_11d.esp5_13.polfac.branchfac<-gof(ev.4_11d.esp5_13.polfac.branchfac)
plot(evgof.4_11d.esp5_13.polfac.branchfac); evgof.4_11d.esp5_13.polfac.branchfac

ev.4_11d.esp5_13.polfac.branchfac.core5 <- ergm(EVNETN~edges + degree(4:11) + esp(5:13) + nodefactor("political") + nodefactor("branch") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(parallel=12, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.esp5_13.polfac.branchfac.core5)
evgof.4_11d.esp5_13.polfac.branchfac.core5<-gof(ev.4_11d.esp5_13.polfac.branchfac.core5)
plot(evgof.4_11d.esp5_13.polfac.branchfac.core5); evgof.4_11d.esp5_13.polfac.branchfac.core5


# We also wanted to test closure and brokerage tendencies using the ergm.graphlets package. While brokerage can be captured by some of our earlier tests (such as k-star and alt-kstar), Yaveroglu's ergm.graphlet package allows to measure the impact of other structural subgraphs on the probability of tie formation in the network. These terms didn't really proffer any meaningful results either. 

grorbetCov(branch, 6)
grorbitCov(political, 6)
grorbetCov(branch, 9:11)
grorbetCov(political, 9:11)

# Open closure didn't converge even after upping the maxit to 50 

ev.gwesp <- ergm(EVNETN~edges + gwesp(), control=control.ergm(rivanna_cores, MCMLE.maxit = 50, parallel.type="PSOCK"))
summary(ev.gwesp) 

# Examining open closure and brokerage alongside current best model 

ev.4_11d.gwesp.kpath <- ergm(EVNETN~edges + degree(4:11) + gwesp() + kpath(), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.gwesp.kpath) 

ev.4_11d.bkbranch6 <- ergm(EVNETN~edges + degree(4:11) + grorbitFactor("branch", 6, 0), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.bkbranch6)
evgof.4_11d.bkbranch6<-gof(ev.4_11d.bkbranch6)
plot(evgof.4_11d.bkbranch6); evgof.4_11d.bkbranch6

ev.4_11d.bkbranch9_11 <- ergm(EVNETN~edges + degree(4:11) + grorbitFactor("branch", 9:11, 0), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.bkbranch9_11)
evgof.4_11d.bkbranch9_11 <- gof(ev.4_11d.bkbranch9_11)
plot(evgof.4_11d.bkbranch9_11); evgof.4_11d.bkbranch9_11

ev.4_11d.bkparty6 <- ergm(EVNETN~edges + degree(4:11) + grorbitFactor("political", 6), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.bkparty6)
evgof.4_11d.bkparty6<-gof(ev.4_11d.bkparty6)
plot(evgof.4_11d.bkparty6); evgof.4_11d.bkparty6

ev.4_11d.bkparty_pol_9_11 <- ergm(EVNETN~edges + degree(4:11) + grorbitFactor("political", 9:11, 0), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.bkparty_pol_9_11)
evgof.4_11d.bkparty_pol_9_11<-gof(ev.4_11d.bkparty_pol_9_11)
plot(evgof.4_11d.bkparty_pol_9_11); evgof.4_11d.bkparty_pol_9_11

# As the result of none of those tests offering anything that helped our model fit, we decided to stick to a structural model predicated on degree distribution without the additional parameters that are useful when nodal attributes are removed completely. This was the best fitting model: 

ev.4_11d <- ergm(EVNETN~edges + degree(4:11), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d)
evgof.4_11d<-gof(ev.4_11d)
plot(evgof.4_11d); evgof.4_11d

# From there, I honed the final models based on politcal party, political branch and subcommittee membership. While both gender and political tenure were included in earlier models, those parameters had very little impact on improving the final models. Thus, we decided to exclude those factors to be most parsimonious. 

ev.4_11d.covten.aff.pol.core5 <- ergm(EVNETN~edges + degree(4:11) + nodecov("tenure") + nodematch("affil") + nodematch("political") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.covten.aff.pol.core5)
evgof.4_11d.covten.aff.pol.core5<-gof(ev.4_11d.covten.aff.pol.core5)
plot(evgof.4_11d.covten.aff.pol.core5); evgof.4_11d.covten.aff.pol.core5

ev.4_11d.aff.core5 <- ergm(EVNETN~edges + degree(4:11) + nodematch("affil") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.aff.core5)
evgof.4_11d.aff.core5<-gof(ev.4_11d.aff.core5)
plot(evgof.4_11d.aff.core5); evgof.4_11d.aff.core5

ev.4_11d.pol.pnp.core5 <- ergm(EVNETN~edges + degree(4:11) + nodematch("political") + nodematch("polnotpol") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.pol.pnp.core5)
EVNETNgof.4_11d.pol.pnp.core5<-gof(ev.4_11d.pol.pnp.core5)
plot(EVNETNgof.4_11d.pol.pnp.core5); EVNETNgof.4_11d.pol.pnp.core5

ev.pol.pnp.core5 <- ergm(EVNETN~edges +  nodematch("political") + nodematch("polnotpol") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.pol.pnp.core5)
evgof.pol.pnp.core5<-gof(ev.pol.pnp.core5)
plot(evgof.pol.pnp.core5); evgof.pol.pnp.core5

ev.4_11d.pol.pnp <- ergm(EVNETN~edges + degree(4:11) + nodematch("political") + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.pol.pnp)
evgof.4_11d.pol.pnp<-gof(ev.4_11d.pol.pnp)
plot(evgof.4_11d.pol.pnp); evgof.4_11d.pol.pnp

ev.4_11d.polfac.pnp <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodematch("polnotpol"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.pnp)
evgof.4_11d.polfac.pnp<-gof(ev.4_11d.polfac.pnp)
plot(evgof.4_11d.polfac.pnp); evgof.4_11d.polfac.pnp

ev.4_11d.polfac.pnp.core5 <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodematch("polnotpol") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.pnp.core5)
evgof.4_11d.polfac.pnp.core5<-gof(ev.4_11d.polfac.pnp.core5)
plot(evgof.4_11d.polfac.pnp.core5); evgof.4_11d.polfac.pnp.core5

ev.polfac.branchfac.core5 <- ergm(EVNETN~edges + nodefactor("political") + nodefactor("branch") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.polfac.branchfac.core5)
evgof.polfac.branchfac.core5<-gof(ev.polfac.branchfac.core5)
plot(evgof.polfac.branchfac.core5); evgof.polfac.branchfac.core5

ev.4_11d.polfac.branchfac <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodefactor("branch"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.branchfac)
evgof.4_11d.polfac.branchfac<-gof(ev.4_11d.polfac.branchfac)
plot(ev.4_11d.polfac.branchfac); ev.4_11d.polfac.branchfac

ev.4_11d.polfac.branchfac.core5 <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodefactor("branch") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.branchfac.core5)
evgof.4_11d.polfac.branchfac.core5<-gof(ev.4_11d.polfac.branchfac.core5)
plot(evgof.4_11d.polfac.branchfac.core5); evgof.4_11d.polfac.branchfac.core5

ev.4_11d.polfac.affilfac.core5 <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodefactor("affil") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.affilfac.core5)
evgof.4_11d.polfac.affilfac.core5<-gof(ev.4_11d.polfac.affilfac.core5)
plot(evgof.4_11d.polfac.affilfac.core5); evgof.4_11d.polfac.affilfac.core5

#########################################################################
# Section 5: Summary of final models   
#########################################################################

# In the end, we chose to include these four models in the final manuscript: 

# Final model 1: Looks at Structural Sweetspot of 4:11 degrees of connectivity 

ev.4_11d <- ergm(EVNETN~edges + degree(4:11), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d)

# Final model 2: Looks at best fitting model of nodal attributes 

ev.polfac.branchfac.core5 <- ergm(EVNETN~edges + nodefactor("political") + nodefactor("branch") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.polfac.branchfac.core5)

# Final model 3: Looks at structural sweetspot in combination with political party and House/Senate membership  

ev.4_11d.polfac.branchfac <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodefactor("branch"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.branchfac)

# Final model 4: Looks at structural sweetspot in combination with political party, branch membership, and impact of top-5 committees 

ev.4_11d.polfac.branchfac.core5 <- ergm(EVNETN~edges + degree(4:11) + nodefactor("political") + nodefactor("branch") + nodematch("appropriations") + nodematch("help") + nodematch("agriculture") + nodematch("energycomm") + nodematch("educ_wf"), control=control.ergm(rivanna_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(ev.4_11d.polfac.branchfac.core5)

# Final model summaries and GOF statistics 

final_m1 <- ev.4_11d
final_m2 <- ev.polfac.branchfac.core5
final_m3 <- ev.4_11d.polfac.branchfac
final_m4 <- ev.4_11d.polfac.branchfac.core5
best_overfit <- ev.4_6d.8_12d.20d.esp5_1319

summary(final_m1)
summary(final_m2)
summary(final_m3)
summary(final_m4)
summary(best_overfit)

gof_m1 <- gof(final_m1) 
gof_m2 <- gof(final_m2) 
gof_m3 <- gof(final_m3) 
gof_m4 <- gof(final_m4) 
gof_overfit <- gof(best_overfit)

plot(gof_m1); gof_m1
plot(gof_m2); gof_m2
plot(gof_m3); gof_m3
plot(gof_m4); gof_m4
plot(gof_overfit); gof_overfit