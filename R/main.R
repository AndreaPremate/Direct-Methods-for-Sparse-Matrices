# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca

# Packages ----------------------------------------------------------------------

# Packages documentation:
# https://cran.r-project.org/web/packages/here/index.html
# https://cran.r-project.org/web/packages/Matrix/index.html
# https://cran.r-project.org/web/packages/bench/index.html

suppressMessages({
  chooseCRANmirror(ind = 1) # select "Cloud" mirror
  packages <- c("here", "bench", "Matrix", "SparseM") # packages to be installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  invisible(lapply(packages, library, character.only = TRUE, quietly = T))
})
rm(list = c("packages", "installed_packages")) # clean memory
options(warn = -1) # clean output

#Initializations ---------------------------------------------------------------

# Variables initialization
matrices_dir <- "matrices/final/" # set matrices directory
list_matrices_mtx <- list.files(path = matrices_dir, pattern = ".mtx$") # get list of .mtx files in the directory
list_matrices <- sub(".mtx$", "", list_matrices_mtx) # get matrices names (for loop purpose)
results_read_csv <- file.path("R/results/read.csv") # set results directory and file for read
results_solve_csv <- file.path("R/results/solve.csv") # set results directory and file for solve
mtx_read_1row <- data.frame("Matrix Name", "Loading Time (s)", "Size (B)", "Memory Allocation (B)")
mtx_solve_1row <- data.frame("Matrix Name", "Total Execution Time (s)", "Relative Error (norm2)", "Total RAM (B)", "Peak RAM (B)")

# CSV Initialization
invisible(file.remove(results_read_csv)) # reset csv files on every start
invisible(file.remove(results_solve_csv)) # reset csv files on every start
close(file(results_read_csv)) # close connection
close(file(results_solve_csv)) # close connection
write.table(mtx_read_1row, results_read_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mtx_solve_1row, results_solve_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
rm(list = c("mtx_read_1row", "mtx_solve_1row")) # clean memory
invisible(gc()) # clean memory
close(file(results_read_csv)) # close connection
close(file(results_solve_csv)) # close connection

# Matrix Load and Solve ----------------------------------------------------------

for (i in seq_along(list_matrices)) {
  # Matrix Read
  cat("\n==================================================\n")
  cat("Reading", list_matrices[i], "matrix:")
  matrix_name <- sub(".mtx$", "", list_matrices_mtx[i]) # get current matrix name
  matrix_dir <- paste0(matrices_dir, list_matrices_mtx[i]) # get current matrix path
  matrix_read <- mark(readMM(matrix_dir), time_unit = "s") # read matrix

  # Matrix Read Results
  assign(paste0(matrix_name), matrix_read$result[[1]]) # assign name to matrix (same as the filename)
  matrix_size <- as.numeric(object.size(matrix_read$result[[1]])) # get matrix size
  matrix_allocation <- as.numeric(matrix_read$mem_alloc) # get matrix memory allocation
  matrix_load_time <- as.numeric(matrix_read$total_time) # get matrix load time (s)

  # Matrix Read Results in CSV
  row <- data.frame(list_matrices[i], matrix_load_time, matrix_size, matrix_allocation)
  write.table(row, results_read_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(file(results_read_csv)) # close connection

  # Console Output (Human Readable Format) and Memory Clean
  cat("\n Matrix size (MB):", matrix_size / (2^20))
  cat("\n Reading memory allocation (MB):", matrix_allocation / (2^20))
  cat("\n Loading time (s):", matrix_load_time)
  cat("\nResults written in", results_read_csv)
  rm(list = c("matrix_dir", "matrix_read", "matrix_allocation", "matrix_load_time", "matrix_size", "row")) # clean memory
  invisible(gc()) # clean memory

  cat("\n==================================================")
  cat("\nSolving", list_matrices[i], "matrix:")

  # Matrix Solve Initialization
  A <- get(list_matrices[i]) # get matrix A
  xe <- rep(1, nrow(A)) # calculate vector xe
  b <- A %*% xe # calculate b

  # Matrix Solve
  x <- mark(solve(A, b, sparse=TRUE, tol=.Machine$double.eps), time_unit = "s") # solve Ax=b (Cholesky or LU)

  # Matrix Solve Results
  time_x <- as.numeric(x$total_time) # x execution time (s)
  relative_error <- norm(x$result[[1]] - xe, "2") / norm(xe, "2") # relative error (norm2)
  peak_ram_x <- as.numeric(max(x$memory[[1]]$bytes)) # x peak ram usage (B)
  tot_ram_x <- as.numeric(sum(x$memory[[1]]$bytes)) # x total ram usage (B)

  # Matrix Solve Results in CSV
  row <- data.frame(list_matrices[i], time_x, relative_error, tot_ram_x, peak_ram_x)
  write.table(row, results_solve_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(file(results_solve_csv)) # close connection

  # Console Output (Human Readable Format) and Memory Clean
  cat("\n Total Execution Time (min):", time_x / 60)
  cat("\n Peak RAM (MB):", peak_ram_x / (2^20))
  cat("\n Relative Error (norm2):", relative_error)
  cat("\nResults written in", results_solve_csv)
  rm(list = c("x", "A", "R", "xe", "b", "matrix_name", "time_x", "relative_error", "tot_ram_x", "peak_ram_x", "row", list_matrices[i]))
  invisible(gc())
}
cat("\n==================================================")
cat("\nMatrices solving finished!")