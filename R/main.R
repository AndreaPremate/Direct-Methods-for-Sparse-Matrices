# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
#package names
packages <- c("here","stringr","pryr","peakRAM","parallel","doSNOW","Matrix","spam","spam64")
options(warn=-1)
# install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# load packages
invisible(lapply(packages, library, character.only = TRUE))

# Reading Matrices (mtx) multicore version ---------------------------------------
cat("\n==================================================\n")

matrices_dir <- "matrices/positive/test_matrix/" #matrices directory
list_matrices_mtx <- list.files(path=matrices_dir, pattern=".mtx$") #list of .mtx files in the directory
list_matrices <- sub(".mtx$", "", list_matrices_mtx) #matrices names for loop purpose
num_matrices <- length(list_matrices) #number of matrices for loop purpose

# Matrix reading loop (automatic name generation and matrix assignation)
for (i in 1:num_matrices){
    matrix_name <- sub(".mtx$", "", list_matrices_mtx[i])
    matrix_dir <- paste0(matrices_dir, list_matrices_mtx[i], sep="")
    matrix_read <- readMM(matrix_dir)
    assign(paste0(matrix_name), matrix_read)
    cat(matrix_name, "memory size: ")
    print(object_size(matrix_read))
}

# Data preparation for loops
results <- data.frame(df_matrix=character(), df_time=double(), df_relative_error=double(), df_ram=double())

# Solving Matrices ---------------------------------------------------------------

for(i in 1:num_matrices){
  cat("\n\n==================================================\n")
  cat("Solving", list_matrices[i], "matrix :")

  # compute B (such that the exact solution of Ax=b is xe=[1,1,..])
  A <- get(list_matrices[i])
  xe <- rep(1,nrow(A))
  b <- A %*% xe

  # solve Ax=b
  start_time <- proc.time()
  ram_used <- peakRAM(x <- solve(A, b))
  execution_time <- proc.time()-start_time

  # time (sec)
  time <- execution_time[["elapsed"]]
  cat("\n Execution Time (s): ", time)

  # relative error (norm2)
  relative_error <- norm(x-xe, type= "2")/norm(xe, type="2")
  cat("\n Relative Error (norm2): ", relative_error)

  # Peak RAM used (MiB)
  ram <- ram_used$Peak_RAM_Used_MiB
  cat("\n RAM Used (MiB): ", ram)

  # Results in dataframe
  result_matrix <- list(df_matrix=list_matrices[i], df_time=time, df_relative_error=relative_error, df_ram=ram)
  results <- rbind(results, result_matrix, stringsAsFactors = FALSE)
}

# Metrics CSV output -------------------------------------------------------------

results_dir <- "R/results/"
write.csv(results, file.path(results_dir, "results.csv"), row.names = FALSE)