# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
cat("\nInitializing the environment...")
options(tryCatchLog.include.full.call.stack = FALSE)
options(tryCatchLog.include.compact.call.stack = FALSE)
suppressMessages({
  chooseCRANmirror(ind=1)
  #package names
  packages <- c("here","bench","tryCatchLog","Matrix")
  # install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  # load packages
  invisible(lapply(packages, library, character.only = TRUE, quietly=T))
})
# free memory
rm(list=c("packages", "installed_packages"))
invisible(gc())

# Variables initialization --------------------------------------------------------
matrices_dir <- "matrices/positive/a/" #matrices directory
results_dir <- "R/results/" #results directory
list_matrices_mtx <- list.files(path=matrices_dir, pattern=".mtx$") #list of .mtx files in the directory
list_matrices <- sub(".mtx$", "", list_matrices_mtx) #matrices names for loop purpose
non_sym_matrices <- vector()

# CSV Initialization --------------------------------------------------------------
mtx_read_1row <- data.frame("Matrix Name", "Loading Time", "Size", "Memory Allocation")
write.table(mtx_read_1row, file.path(results_dir, "results_matrices_readed.csv"), sep=",", append=TRUE, quote=FALSE, row.names = FALSE, col.names=FALSE)
mtx_solve_1row <- data.frame("Matrix Name", "Solving Time", "Relative Error", "Peak RAM")
write.table(mtx_solve_1row, file.path(results_dir, "results_matrices_solved.csv"), sep=",", append=TRUE, quote=FALSE, row.names = FALSE, col.names=FALSE)

# Matrices PRE-Processing (reading, non-symmetrical matrix deletion) -------------
cat("\nReading the matrices...")

# Matrix reading loop (automatic name generation and matrix assignation)
for (i in seq_along(list_matrices)){
  cat("\n==================================================\n")
  cat("Reading", list_matrices[i], "matrix:")
  matrix_name <- sub(".mtx$", "", list_matrices_mtx[i])
  matrix_dir <- paste0(matrices_dir, list_matrices_mtx[i])
  matrix_read_benchmark <- mark(as(matrix_read <- readMM(matrix_dir), "CsparseMatrix"))

  if (isSymmetric(matrix_read)){
    assign(paste0(matrix_name), matrix_read) #assign matrix to names
    cat("\n Matrix size (Bytes): ")
    matrix_size <- as.numeric(object.size(matrix_read))
    cat(matrix_size)
    cat("\n Reading memory allocation (Bytes): ")
    matrix_allocation <- as.numeric(matrix_read_benchmark$mem_alloc)
    cat(matrix_allocation)
    cat("\n Loading time (s): ")
    matrix_load_time <- matrix_read_benchmark$total_time #time in s
    cat(matrix_load_time)

    # Matrix reading results in CSV
    row <- data.frame(df_matrix=list_matrices[i], df_loading_time=matrix_load_time, df_matrix_size=matrix_size, df_memory_allocation=matrix_allocation)
    write.table(row, file.path(results_dir, "results_matrices_readed.csv"), sep=",", append=TRUE, quote=FALSE, row.names = FALSE, col.names=FALSE)
    cat("\nResults written in", paste0(results_dir,"results_matrices_readed.csv"))
  } else {
    cat("\n Matrix non-symmetrical!") #save index of non-symmetrical matrices
    non_sym_matrices[length(non_sym_matrices)+1] <- i
  }
}

cat("\n==================================================")
cat("\nDeleting non-symmetrical matrices... ")
non_sym_matrices <- non_sym_matrices[!is.na(non_sym_matrices)]
list_matrices <- list_matrices[-non_sym_matrices]
cat("\nNon-symmetrical matrices deleted!")
# free memory
rm(list=c("list_matrices_mtx", "i", "matrix_read","matrix_name","matrix_dir", "matrix_read_benchmark","matrix_size","matrix_allocation","matrix_load_time"))
invisible(gc()) #called after rm as suggested in R documentation
cat("\n==================================================")
cat("\nMatrices reading finished!")

# Matrices Processing (solving) --------------------------------------------------
cat("\n\nSolving the matrices...")
for(i in seq_along(list_matrices)){
  cat("\n==================================================\n")
  cat("Solving", list_matrices[i], "matrix:\n")

  # cholesky analysis
  A <- get(list_matrices[i])
  positive_matrix_check <- tryCatchLog(
    { R <- chol(A) },
    error=function(e){ cat("Matrix not-positive definite, deleted!") }
  )

  # solve positive matrix
  if(inherits(positive_matrix_check, "error")){
    xe <- rep(1,nrow(A))
    b <- crossprod(A, xe)

    # solve Ax=b
    x_benchmark <- mark(x <- solve(R, b))

    # time (sec)
    time <- x_benchmark$total_time #time in s
    conv_time <- time/60 #time in min
    cat(" Execution Time (min): ", conv_time)
    # relative error (norm2)
    relative_error <- norm(x-xe,"2")/norm(xe,"2")
    cat("\n Relative Error (norm2): ", relative_error)
    # RAM used (Bytes)
    conv_ram <- as.numeric(x_benchmark$mem_alloc)
    cat("\n RAM used (Bytes): ", conv_ram)

    # Matrix solving results in dataframe
    row1 <- data.frame(df_matrix=list_matrices[i], df_time=conv_time, df_relative_error=relative_error, df_peak_ram=conv_ram)
    write.table(row1, file.path(results_dir, "results_matrices_solved.csv"), sep=",", append=TRUE, quote=FALSE, row.names = FALSE, col.names=FALSE)
    cat("\nResults written in", paste0(results_dir,"results_matrices_solved.csv"))
  }

  # free memory (dealing with big matrix could generate memory errors)
  suppressWarnings(
     rm(list=c("A","xe","b","x_benchmark","time", "conv_time", "relative_error","conv_ram",list_matrices[i])))
  invisible(gc()) #called after rm as suggested in R documentation
}
cat("\n==================================================")
cat("\nMatrices solving finished!")