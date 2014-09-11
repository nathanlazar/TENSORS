#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to run Bayesian Multiview Tensor Factorization package
# on cell line, drug, dosage data from the paper 
# "Modeling precision treatment of breast cancer" by Daemen et al.
# published in Genome Biology 2013

# Usage:
# BC_bmtf_1.0.R
#   <bmtf directory>
#   <working directory>
#   <min K value>
#   <max K value>

# Example:
# BC_CellLines1.0.R
#   /mnt/lustre1/users/lazar/TENSORS/BMTF/bmtf/
#   /mnt/lustre1/users/lazar/TENSORS/BC_CellLines/
#   3
#   30

args <- commandArgs(TRUE)
bmtfdir <- args[1]
outdir <- args[2]
minK <- args[3]
naxK <- args[4]

library(dplyr)
require(compiler)
require(tensor)

# source code files
source(paste0(bmtfdir, "bmtf.R"))
source(paste0(bmtfdir, "tools_data.R"))
source(paste0(bmtfdir, "tools_plot.R"))
source(paste0(bmtfdir, "getRandomInits.R"))
source(paste0(bmtfdir, "tools_test.R"))
source(paste0(outdir, "tensor_functions.R"))

response <- read.table(paste0(outdir, 'gb-2013-14-10-r110-S9.txt'), 
                       header=T)
resp <- tbl_df(response)

# Build the tensor
##############################

# Normalize by baseline od values?
#mutate(resp, )

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

opts <- getDefaultOpts(tens)

# Save data so that workers can access it
save.image(paste0(outdir, 'BC_bmtf.dat'))

# Bayesian mult-view tensor factorization
############################################

# submit HTCondor jobs 
system('condor_submit bmtf_submit')

# Wait until all have finished running
while(!exists(

# Read in output files

for(k in minK:maxK) {

}

K <- 4
model <- bmtf(list(tens),K,opts) # Run the algorithm
model.mean <- getPosteriorMean(model)$rs.mean





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


