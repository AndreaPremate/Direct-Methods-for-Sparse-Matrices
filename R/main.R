# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
#package names
packages <- c("here","stringr","pryr","peakRAM","Matrix","SparseM")
options(warn=-1)
# install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# load packages
invisible(lapply(packages, library, character.only = TRUE))

# Tables initialization for CSV --------------------------------------------------
matrix_reading_results <- data.frame(df_matrix=character(), df_loading_time=double(), df_memory_size=double())
matrix_solving_results <- data.frame(df_matrix=character(), df_time=double(), df_relative_error=double(), df_ram=double())

# Reading Matrices (mtx) ---------------------------------------------------------
matrices_dir <- "matrices/positive/" #matrices directory
list_matrices_mtx <- list.files(path=matrices_dir, pattern=".mtx$") #list of .mtx files in the directory
list_matrices <- sub(".mtx$", "", list_matrices_mtx) #matrices names for loop purpose
num_matrices <- length(list_matrices) #number of matrices for loop purpose

# Matrix reading loop (automatic name generation and matrix assignation)
cat("\n\nReading matrices...\n")
for (i in 1:num_matrices){
    cat("\n==================================================\n")
    cat("Reading", list_matrices[i], "matrix:")
    matrix_name <- sub(".mtx$", "", list_matrices_mtx[i])
    matrix_dir <- paste0(matrices_dir, list_matrices_mtx[i], sep="")
    start_loading_time <- proc.time()
    matrix_read <- readMM(matrix_dir)
    matrix_loading_time <- proc.time()-start_loading_time
    assign(paste0(matrix_name), matrix_read)

    cat("\n Memory size (Bytes): ")
    cat(object_size(matrix_read))
    cat("\n Loading time (s): ")
    cat(matrix_loading_time[["elapsed"]])

    # Matrix reading results in dataframe
    matrix_reading_data <- list(df_matrix=list_matrices[i], df_loading_time=matrix_loading_time[["elapsed"]], df_memory_size=object_size(matrix_read))
    matrix_reading_results <- rbind(matrix_reading_results, matrix_reading_data, stringsAsFactors = FALSE)
}
cat("\n==================================================")
cat("\n\nMatrices reading finished!")

# Solving Matrices ---------------------------------------------------------------
cat("\n\nSolving matrices...\n")
for(i in 1:num_matrices){
  cat("\n==================================================\n")
  cat("Solving", list_matrices[i], "matrix:")

  # compute B (such that the exact solution of Ax=b is xe=[1,1,..])
  A <- get(list_matrices[i])
  xe <- rep(1,nrow(A))
  b <- A %*% xe

  # solve Ax=b
  start_time <- proc.time()
  ram_used <- peakRAM(x <- SparseM::solve(A, b))
  execution_time <- proc.time()-start_time

  # time (sec)
  time <- execution_time[["elapsed"]]
  cat("\n Execution Time (s): ", time)

  # relative error (norm2)
  relative_error <- norm(x-xe,"2")/norm(xe,"2")
  cat("\n Relative Error (norm2): ", relative_error)

  # Peak RAM used (MiB)
  ram <- ram_used$Peak_RAM_Used_MiB
  cat("\n RAM Used (MiB): ", ram)

  # Matrix solving results in dataframe
  matrix_solving_data <- list(df_matrix=list_matrices[i], df_time=time, df_relative_error=relative_error, df_ram=ram)
  matrix_solving_results <- rbind(matrix_solving_results, matrix_solving_data, stringsAsFactors = FALSE)
}
cat("\n==================================================")
cat("\n\nMatrices solving finished!")

# Metrics CSV output -------------------------------------------------------------
results_dir <- "R/results/"
write.csv(matrix_reading_results, file.path(results_dir, "matrix_reading_results.csv"), row.names = FALSE)
write.csv(matrix_solving_results, file.path(results_dir, "matrix_solving_results.csv"), row.names = FALSE)
cat("\n\nResults CSV file written!")