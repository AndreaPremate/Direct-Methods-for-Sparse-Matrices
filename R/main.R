# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
install.packages("pacman")
library(pacman)
pacman::p_load(here, parallel, MASS, Matrix)

# Reading Matrices (mtx) ---------------------------------------------------------
hook_1498 <- readMM('matrices/positive/Hook_1498.mtx')
g3_circuit <- readMM('matrices/positive/G3_circuit.mtx')
nd24k <- readMM('matrices/positive/nd24k.mtx')
bundle_adj <- readMM('matrices/positive/bundle_adj.mtx')
ifiss_mat <- readMM('matrices/non_positive/ifiss_mat.mtx')
tsc_opf_1047 <- readMM('matrices/non_positive/tsc_opf_1047.mtx')
ns3Da <- readMM('matrices/non_positive/ns3Da.mtx')
gt01r <- readMM('matrices/non_positive/gt01r.mtx')
# matrices are already in sparse format so they don't need to be converted

# Memory (RAM) occupied by each matrix
cat("\n==================================================\n")
cat("Hook_1498 memory usage: ")
print(object.size(hook_1498), units="Mb")
cat("G3_circuit memory usage: ")
print(object.size(g3_circuit), units="Mb")
cat("nd24k memory usage: ")
print(object.size(nd24k), units="Mb")
cat("bundle_adj memory usage: ")
print(object.size(bundle_adj), units="Mb")
cat("ifiss_mat memory usage: ")
print(object.size(ifiss_mat), units="Mb")
cat("TSC_OPF_1047 memory usage: ")
print(object.size(tsc_opf_1047), units="Mb")
cat("ns3Da memory usage: ")
print(object.size(ns3Da), units="Mb")
cat("GT01R memory usage: ")
print(object.size(gt01r), units="Mb")

# Solving Matrices ---------------------------------------------------------------
# Compute B (such that the exact solution of Ax=b is xe=[1,1,..])
a <- gt01r
xe <- rep(1,nrow(a))
b <- a*xe

# Solve Ax=b
x <- solve(a,b)