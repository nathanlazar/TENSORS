#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to run Bayesian Multiview Tensor Factorization package
# on cell line, drug, dosage data from the paper 
# "Modeling precision treatment of breast cancer" by Daemen et al.
# published in Genome Biology 2013

# Usage:
# BC_CellLines1.0.R
#   <bmtf directory>
#   <working directory>
#   <min K value>
#   <max K value>

# Example:
# BC_CellLines1.0.R
#   /mnt/lustre1/users/lazar/TENSORS/BMTF/
#   /mnt/lustre1/users/lazar/TENSORS/BC_CellLines/
#   3
#   30

args <- commandArgs(TRUE)
bindir <- args[1]
outdir <- args[2]
minK <- args[3]
naxK <- args[4]

library(xlsx)
library(dplyr)
require(compiler)
require(tensor)
# Read in the code
setwd('C:/Users/lazar/Desktop/Thesis/Tensors/BMTF/bmtf/')
source("tools_test.R")
source("bmtf.R")
setwd("C:/Users/lazar/Desktop/Thesis/BC_CellLines")

response <- read.table('data/gb-2013-14-10-r110-S9.txt', header=T)
resp <- tbl_df(response)
gi50 <- read.xlsx('data/gb-2013-14-10-r110-S4.xlsx', 1,
                  startRow=3)

plot_curves <- function(resp, thresh, n) {
  dose <- unlist(lapply(0:9, rep, 3))
  for(i in 1:n) {
    plot(dose, select(resp[i,], od0.1:od9.3),
         main=paste('Cell Line:',resp[i,1], 
                    'Drug:', resp[i,2]), 
         xlab="concentration",
         ylab='cells')
    lines(c(thresh, thresh), c(0,10000))
  }
}

#plot the drug response data with GI50 cutoffs indicated
for(i in 1:nrow(gi50)) {
  drug <- as.character(gi50[i,1])
  thresh <- gi50[i,3]
  sub <- filter(resp, compound==drug)
#  print(sub)
  if(nrow(sub)>0) {
    plot_curves(sub, thresh, 5)
  }
}

# Build the tensor
##############################

# Normalize by baseline od values?
#mutate(resp, )

# Function to use dplyr to average responses for each dose
avg_resp <- function(resp, cl, cpd) {
tmp <- resp %.% 
  filter(cellline==cl) %.%
  filter(compound==cpd) %.%
  select(od0.1:od9.3) 
  if(nrow(tmp)>0) {
    tmp <- tmp %.%
    summarise_each(funs(mean)) %.%
    mutate(od0 = mean(od0.1, od0.2, od0.3),
           od1 = mean(od1.1, od1.2, od1.3),
           od2 = mean(od2.1, od2.2, od2.3),
           od3 = mean(od3.1, od3.2, od3.3),
           od4 = mean(od4.1, od4.2, od4.3),
           od5 = mean(od5.1, od5.2, od5.3),
           od6 = mean(od6.1, od6.2, od6.3),
           od7 = mean(od7.1, od7.2, od7.3),
           od8 = mean(od8.1, od8.2, od8.3),
           od9 = mean(od9.1, od9.2, od9.3)) %.%
    select(od0:od9)
    return(as.numeric(tmp))
  } else {
    return(rep(NA, 10))
  }
}

# Arrange data into a tensor 
# (rows are cell lines, columns are drugs, 3rd dim are dosages)
# There are 71 cell lines, 90 drugs and 10 dosages
tens <- array(NA, c(71,90,10))
dimnames(tens) <- list(levels(resp$cellline),
                       levels(resp$compound), 
                       paste0('mean.od', 0:9))
for(cellline in levels(resp$cellline)) {
  for(drug in levels(resp$compound))
    tens[cellline, drug,] <- avg_resp(resp, cellline, drug)
}

# Bayesian mult-view tensor factorization
############################################

K <- 4
opts <- getDefaultOpts(tens)
model <- bmtf(list(tens),K,opts) # Run the algorithm
model.mean = getPosteriorMean(model)$rs.mean





###########################################################
##
## Generate some data from the model, with pre-specified
## latent components
##
M = 2
N = 30
D = c(40,50)
L = 20
SNR = c(3,2)
CP <- getVisualComponents.3K.Paper(N,N,L,D,M)
tensor.data <- data.train.test.1(CP,SNR) 

K <- ncol(tensor.data$Z) + 1 # K is the number of components
opts <- getDefaultOpts(tensor.data$Y)
model <- bmtf(tensor.data$Y,K,opts) # Run the algorithm

tensor.data
length(tensor.data) #length 17
length(tensor.data[[3]])
dim(tensor.data)


