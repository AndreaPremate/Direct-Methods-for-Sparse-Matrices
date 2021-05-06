# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
#package names
packages <- c("here","stringr","pryr","peakRAM","parallel", "doSNOW", "Matrix")
options(warn=-1)
# install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# load packages
invisible(lapply(packages, library, character.only = TRUE))

# Reading Matrices (mtx) ---------------------------------------------------------
cat("\n==================================================\n")

load_matrices <- function(){
  cl <- parallel::makeCluster(num_cores, setup_strategy = "sequential")
  doParallel:: registerDoParallel(cl)

  hook_1498 <- readMM('matrices/positive/Hook_1498.mtx')
  g3_circuit <- readMM('matrices/positive/G3_circuit.mtx')
  nd24k <- readMM('matrices/positive/nd24k.mtx')
  bundle_adj <- readMM('matrices/positive/bundle_adj.mtx')
  ifiss_mat <- readMM('matrices/non_positive/ifiss_mat.mtx')
  tsc_opf_1047 <- readMM('matrices/non_positive/tsc_opf_1047.mtx')
  ns3Da <- readMM('matrices/non_positive/ns3Da.mtx')
  gt01r <- readMM('matrices/non_positive/gt01r.mtx')
  # matrices are already in sparse format so they don't need to be converted

  parallel::stopCluster(cl)
}
load_matrices()

# data preparation for loops
matrices <- c("hook_1498", "g3_circuit", "nd24k", "bundle_adj", "ifiss_mat", "tsc_opf_1047", "ns3Da", "gt01r")
n_matrices <- length(matrices)
num_cores <- detectCores()-1

# Memory size of each matrix (multicore mode)
#for(i in 1:n_matrices){
#  cat(matrices[i], "memory size: ")
#  matrix<-get(matrices[i])
#  print(object_size(matrix))
#}

#read_matrices <- function(input_matrices, i){
#  #cat(input_matrices[i], "memory size: ")
#  matrix<-get(input_matrices[i])
#  print(object_size(matrix))
#}
#mclapply(n_matrices, read_matrices(matrices, i), mc.cores = num_cores)

read_matrices <- function(input_matrices, num_matrices){
  cl <- parallel::makeCluster(num_cores, setup_strategy = "sequential")
  doParallel:: registerDoParallel(cl)

  for(i in 1:num_matrices){
  cat(input_matrices[i], "memory size: ")
  matrix<-get(input_matrices[i])
  print(object_size(matrix))
  }

  parallel::stopCluster(cl)
}
read_matrices(matrices, n_matrices)


# Solving Matrices ---------------------------------------------------------------

# compute B (such that the exact solution of Ax=b is xe=[1,1,..])
A <- gt01r
xe <- rep(1,nrow(A))
b <- A*xe

# solve Ax=b
start_time <- proc.time()
ram_used <- peakRAM(x <- solve(A, b))
execution_time <- proc.time()-start_time

# Metrics ------------------------------------------------------------------------
cat("\n==================================================\n")

# time (rbenchmark/microbenchmart not used because executions can last over 5 minute)
cat("Execution Time (sec): ")
cat(execution_time[["elapsed"]], "\n")
# elapsed time = CPU execution time + external processes (OS)

# relative error
#cat("Relative Error (norm2): ")
#relative_error <- norm(x-xe, "2")/norm(xe, "2")
#cat(relative_error, "\n")

# RAM used
cat("RAM used (MiB): ")
cat(ram_used$Peak_RAM_Used_MiB, "\n")
# I do not extract the elapsed time of this function because it does not consider the time that elapses between the external processes of the OS (it considers only CPU execution time)


