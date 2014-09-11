#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to run Bayesian Multiview Tensor Factorization package
# in parallel for different values of K given data saved in 
# BC_bmtf.dat file

# Usage: ./run_bmtf.R
#   /mnt/lustre1/users/lazar/TENSORS/BC_CellLines/BC_bmtf.dat
#   5

require(compiler)
require(tensor)

args <- commandArgs(TRUE)
datafile <- args[1]
tmp <- strsplit(datafile, '/')[[1]]
wdir <- paste(tmp[1:(length(tmp)-1)], collapse='/')
wdir <- paste0(wdir, '/')
K <- as.numeric(args[2])

# Load the file saved by BC_bmtf_1.0.R which should have the tensor 
# 'tens', an options object 'opts'  as well as the functions needed 
# to run the BMTF model
load(datafile)

model <- bmtf(list(tens), K, opts)
model.mean <- getPosteriorMean(model)$rs.mean

save(model.mean, file=paste0(wdir, 'model_', K, '.dat'))