#############################################################################3
# SAMPLE SIZE IN METAGENOMICS: DETERMINATION OF
# RICHNESS AND SAMPLE SIZE (EFFORT)
# USING A BAYESIAN APPROACH
##############################################################################
# Toni Monleon-Getino, Clara I Rodriguez-Casado (Section of Statistics, University of Barcelona).
# Group of Research in Biostatistics and Bioinformatics
# amonleong@ub.edu
# library Biost3. 30-8-2017



#see examples of JAGS (Bayes MCMC) at: http://rstudio-pubs-static.s3.amazonaws.com/15236_9bc0cd0966924b139c5162d7d61a2436.html


##################################################################
# EXAMPLES OF USE OF THE FUNCTION MetagenSample.size.H1()
# USING SIMULATION OF A METAGENOMIC DISTRIBUTION
# CALCULATION OF THE RICHNESS-EFFORT OF A ABUNDANCE-RICHNESS METAGENOMIC MATRIX (ROWS:OTUs, COLUMNS:SAMPLES, CELSS: Abundance) )
###################################################################

###################################################
# 1) SIMULATE A OVER-SAMPLING METAGENOMIC MATRIX
###################################################

# PARAMETERS FOR THE SIMULATION
nsites<-100 #NUMBER OF OTUs
nsimulac<-15 #SITES OR SAMPLES (REPLICATIONS)
ab.total<-300 #ABUNDANCE IN EACH SITE
library(LearnBayes)
ppp <- rdirichlet(1, par = rep(1, nsites))
ppp
X
N
matriu <- array(0, dim=c(nsites, nsimulac))
for(i in 1:nsimulac){
  X <- as.vector(rmultinom(1, size =ab.total , prob = ppp))
  matriu[,i]<-X
  N <- sum(X)
}
#SIMULATED METAGENOMIC-MATRIX
matriu #


#COMPUTE THE RICHNESS AND SAMPLING EFFORT BAYESIAN APPROACH
#Upload the library
library(BDSbiost3)
MetagenSample.size.H1(matriu, quart.cut=0.3,type=T,model.probability=1) #ejemplo de no saturacion

#ESTIMATION OF THE ASYMPTOTIQUE RICHNESS USING THE FUNCTION
# poolaccum package Vegan for R
## Accumulation model
library(vegan)
plot(poolaccum(t(matriu))) #plot de Richness
poolaccum(t(matriu))




###################################################
# 2) SIMULATE A UNDER-SAMPLING METAGENOMIC MATRIX
###################################################

# PARAMETERS FOR THE SIMULATION
nsites<-100 #NUMBER OF OTUs
nsimulac<-4 #SITES OR SAMPLES (REPLICATIONS)
ab.total<-300 #ABUNDANCE IN EACH SITE
library(LearnBayes)
ppp <- rdirichlet(1, par = rep(1, nsites))
ppp
X
N
matriu <- array(0, dim=c(nsites, nsimulac))
for(i in 1:nsimulac){
  X <- as.vector(rmultinom(1, size =ab.total , prob = ppp))
  matriu[,i]<-X
  N <- sum(X)
}
#SIMULATED METAGENOMIC-MATRIX
matriu #


#COMPUTE THE RICHNESS AND SAMPLING EFFORT BAYESIAN APPROACH
#Upload the library
library(BDSbiost3)
MetagenSample.size.H1(matriu, quart.cut=0.3,type=T,model.probability=1) #ejemplo de no saturacion

#ESTIMATION OF THE ASYMPTOTIQUE RICHNESS USING THE FUNCTION
# poolaccum package Vegan for R
## Accumulation model
library(vegan)
plot(poolaccum(t(matriu))) #plot de Richness
poolaccum(t(matriu))
