# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
install.packages("pacman")
library(pacman)
pacman::p_load(here, Matrix)

# Reading Matrices (mtx) ---------------------------------------------------------
bundleAdj <- readMM('matrices/positive/bundle_adj.mtx')
g3Circuit <- readMM('matrices/positive/G3_circuit.mtx')
hook1498 <- readMM('matrices/positive/Hook_1498.mtx')
nd24k <- readMM('matrices/positive/nd24k.mtx')
# The matrices are already in sparse format so they do not need to be converted.

# Memory (RAM) required from every matrix
object.size(bundleAdj)
object.size(g3Circuit)
object.size(hook1498)
object.size(nd24k)