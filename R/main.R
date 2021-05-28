# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
cat("\nInitializing the environment...")

# Packages documentation:
# https://cran.r-project.org/web/packages/here/index.html
# https://cran.r-project.org/web/packages/Matrix/index.html
# https://cran.r-project.org/web/packages/bench/index.html
#
# Auto-Mirror + Install only packages not yet installed
suppressMessages({
  chooseCRANmirror(ind = 1)
  packages <- c("here", "bench", "Matrix")
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  invisible(lapply(packages, library, character.only = TRUE, quietly = T))
})
rm(list = c("packages", "installed_packages"))

#Initializations ---------------------------------------------------------------
options(warn = -1)

# Variables initialization
matrices_dir <- "matrices/final/" # set matrices directory
list_matrices_mtx <- list.files(path = matrices_dir, pattern = ".mtx$") # get list of .mtx files in the directory
list_matrices <- sub(".mtx$", "", list_matrices_mtx) # get matrices names (for loop purpose)
results_dir <- "R/results/" # set results directory
results_read_csv <- file.path(results_dir, "results_read.csv") # set results directory
results_solve_csv <- file.path(results_dir, "results_solve.csv") # set results directory
mtx_read_1row <- data.frame("Matrix Name", "Loading Time (s)", "Size (Bytes)", "Memory Allocation (Bytes)")
mtx_solve_1row <- data.frame("Matrix Name", "Cholesky Solving Time (s)", "x Solving Time (s)", "Total Solving Time (s)", "Relative Error", "RAM solving Cholesky (Bytes)", "Peak RAM solving Cholesky (Bytes)", "RAM solving x (Bytes)", "Peak RAM solving x (Bytes)")

# CSV Initialization
invisible(file.remove(results_read_csv)) # delete csv files on every start
invisible(file.remove(results_solve_csv)) # delete csv files on every start
write.table(mtx_read_1row, results_read_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mtx_solve_1row, results_solve_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

# Matrix Load and Solve ----------------------------------------------------------
for (i in seq_along(list_matrices)) {
  # Matrix Read
  cat("\n==================================================\n")
  cat("Reading", list_matrices[i], "matrix:")
  matrix_name <- sub(".mtx$", "", list_matrices_mtx[i]) # get current matrix name
  matrix_dir <- paste0(matrices_dir, list_matrices_mtx[i]) # get current matrix path
  matrix_read <- mark(readMM(matrix_dir), time_unit = "s") # read matrix

  # Symmetry Check
  if (isSymmetric(matrix_read$result[[1]])) {
    # Matrix Load
    assign(paste0(matrix_name), matrix_read$result[[1]]) # assign name to matrix (same as the filename)
    matrix_size <- as.numeric(object.size(matrix_read$result[[1]])) # get matrix size
    matrix_allocation <- as.numeric(matrix_read$mem_alloc) # get matrix memory allocation
    matrix_load_time <- as.numeric(matrix_read$total_time) # get matrix load time (s)

    # Results in CSV
    row <- data.frame(list_matrices[i], matrix_load_time, matrix_size, matrix_allocation)
    write.table(row, results_read_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

    # Console Output and Memory Clean
    cat("\n Matrix size (MB):", matrix_size/(2^20))
    cat("\n Reading memory allocation (MB):", matrix_allocation/(2^20))
    cat("\n Loading time (s):", matrix_load_time)
    cat("\nResults written in", paste0(results_dir, "results_read.csv"))
    rm(list = c("matrix_dir", "matrix_read", "matrix_allocation", "matrix_load_time", "row")) # clean memory
    invisible(gc()) # clean memory
  } else {
    cat("\n Matrix non-symmetrical!")
    rm(list = c("matrix_read", "matrix_name", "matrix_dir", "matrix_load_time")) # clean memory
    invisible(gc()) # clean memory
    cat("\n Deleted, skipping to next one...")
    next
  }

  cat("\n==================================================")
  cat("\nSolving", list_matrices[i], "matrix:")

  # Matrix Solve
  A <- get(list_matrices[i]) # get matrix A

  # Not-Positive Definite Check (Cholesky Analysis)
  R <- tryCatch({
    mark(Cholesky(A)) # calculate cholesky
    }
    , error = function(e) {
      cat("\n Matrix not-positive definite!")
    })

  if (!inherits(R, "error") & !is.null(R)) {
    # Solve Ax=b
    xe <- rep(1, nrow(A)) # calculate vector xe
    b <- crossprod(A, xe) # calculate b
    x <- mark(solve(R$result[[1]], b), time_unit = "s") # calculate Ax=b

    # Metrics
    time_chol <- as.numeric(R$total_time) # cholesky solution time (s)
    time_x <- as.numeric(x$total_time) # x solution time (s)
    time_total <- as.numeric(time_x + time_chol) # total solution time (s)
    relative_error <- norm(x$result[[1]] - xe, "2") / norm(xe, "2") # relative error (norm2)
    peak_conv_ram_chol <- as.numeric(max(R$memory[[1]]$bytes)) # cholesky peak ram usage (B)
    conv_ram_chol <- as.numeric(sum(R$memory[[1]]$bytes)) # cholesky total ram usage (B)
    peak_conv_ram_x <- as.numeric(max(x$memory[[1]]$bytes)) # x peak ram usage (B)
    conv_ram_x <- as.numeric(sum(x$memory[[1]]$bytes)) # x total ram usage (B)

    # Results in CSV
    row1 <- data.frame(list_matrices[i], time_chol, time_x, time_total, relative_error, conv_ram_chol, peak_conv_ram_chol, conv_ram_x, peak_conv_ram_x)
    write.table(row1, results_solve_csv, sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat("\nResults written in", paste0(results_dir, "results_solve.csv"))

    # Console Output and Memory Clean
    cat("\n Resolution Time (min):")
    cat("\n  Cholesky (min):", time_chol / 60)
    cat("\n  x (min):", time_x / 60)
    cat("\n  Total (min):", time_total / 60)
    cat("\n RAM used (MB):")
    cat("\n  Peak Cholesky (MB): ", peak_conv_ram_chol/(2^20))
    cat("\n  Total Cholesky (MB): ", conv_ram_chol/(2^20))
    cat("\n  Peak x (MB): ", peak_conv_ram_x/(2^20))
    cat("\n  Total x (MB): ", conv_ram_x/(2^20))
    cat("\n Relative Error (norm2): ", relative_error)
    rm(list = c("x", "A", "R", "xe", "b", list_matrices[i]))
    invisible(gc())
  } else {
    rm(list = c("A", "R", "matrix_read", "matrix_name", "matrix_dir", list_matrices[i]))
    invisible(gc())
    cat("\n Deleted, skipping to next one...")
    next
  }
}
cat("\n==================================================")
cat("\nMatrices solving finished!")