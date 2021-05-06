# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
#package names
packages <- c("here","stringr","pryr","peakRAM","Matrix")
# install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# load packages
invisible(lapply(packages, library, character.only = TRUE))
cat("\n==================================================\n")

# Reading Matrices (mtx) ---------------------------------------------------------
#hook_1498 <- readMM('matrices/positive/Hook_1498.mtx')
#g3_circuit <- readMM('matrices/positive/G3_circuit.mtx')
#nd24k <- readMM('matrices/positive/nd24k.mtx')
#bundle_adj <- readMM('matrices/positive/bundle_adj.mtx')
#ifiss_mat <- readMM('matrices/non_positive/ifiss_mat.mtx')
#tsc_opf_1047 <- readMM('matrices/non_positive/tsc_opf_1047.mtx')
#ns3Da <- readMM('matrices/non_positive/ns3Da.mtx')
gt01r <- readMM('matrices/non_positive/gt01r.mtx')
# matrices are already in sparse format so they don't need to be converted


# Memory size of each matrix
#cat("\n==================================================\n")
#cat("Hook_1498 memory usage: ")
#print(object.size(hook_1498), units="Mb")
#cat("G3_circuit memory usage: ")
#print(object.size(g3_circuit), units="Mb")
#cat("nd24k memory usage: ")
#print(object.size(nd24k), units="Mb")
#cat("bundle_adj memory usage: ")
#print(object.size(bundle_adj), units="Mb")
#cat("ifiss_mat memory usage: ")
#print(object.size(ifiss_mat), units="Mb")
#cat("TSC_OPF_1047 size: ")
#print(object.size(tsc_opf_1047), units="Mb")
#cat("ns3Da memory size: ")
#print(object.size(ns3Da), units="Mb")
cat("GT01R memory size: ")
print(object_size(gt01r))

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
cat("Relative Error (norm2): ")
relative_error <- norm(x-xe, "2")/norm(xe, "2")
cat(relative_error, "\n")

# RAM used
cat("RAM used (MiB): ")
cat(ram_used$Peak_RAM_Used_MiB, "\n")
# I do not extract the elapsed time of this function because it does not consider the time that elapses between the external processes of the OS (it considers only CPU execution time)


