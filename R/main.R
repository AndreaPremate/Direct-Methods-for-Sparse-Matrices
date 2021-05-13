# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
cat("\nInitializing the environment...")
options(warn=-1)
suppressMessages({
  chooseCRANmirror(ind=1)
  #package names
  packages <- c("here","bench","Matrix","SparseM")
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

# Variables initializaion --------------------------------------------------------
cat("\nStarting the program...")
matrices_dir <- "matrices/positive/a/" #matrices directory
results_dir <- "R/results/" #results directory
list_matrices_mtx <- list.files(path=matrices_dir, pattern=".mtx$") #list of .mtx files in the directory
list_matrices <- sub(".mtx$", "", list_matrices_mtx) #matrices names for loop purpose
non_sym_matrices <- vector()
matrix_reading_results <- data.frame(df_matrix=character(), df_loading_time=double(), df_matrix_size=double(), df_memory_allocation=double())
matrix_solving_results <- data.frame(df_matrix=character(), df_time=double(), df_relative_error=double(), df_peak_ram=double())

# Matrices PRE-Processing (reading, non-symmetrical matrix deletion) -------------
cat("\nReading the matrices...")

# Matrix reading loop (automatic name generation and matrix assignation)
for (i in seq_along(list_matrices)){
  cat("\n==================================================\n")
  cat("Reading", list_matrices[i], "matrix:")
  matrix_name <- sub(".mtx$", "", list_matrices_mtx[i])
  matrix_dir <- paste0(matrices_dir, list_matrices_mtx[i])
  matrix_read_benchmark <- mark(matrix_read <- readMM(matrix_dir)) #read.MM on mac doesn't work

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

    # Matrix reading results in dataframe
    matrix_reading_data <- list(df_matrix=list_matrices[i], df_loading_time=matrix_load_time, df_matrix_size=matrix_size, df_memory_allocation=matrix_allocation)
    matrix_reading_results <- rbind(matrix_reading_results, matrix_reading_data, stringsAsFactors = FALSE)
  } else {
    cat("\n Matrix non-symmetrical!") #save index of non-symmetrical matrices
    non_sym_matrices[length(non_sym_matrices)+1] <- i
  }

  # free memory (dealing with big matrix could generate memory errors)
  rm(list=c("matrix_read","matrix_name","matrix_dir", "matrix_read_benchmark","matrix_size","matrix_allocation","matrix_load_time","matrix_reading_data"))
  invisible(gc()) #called after rm as suggested in R documentation
}

cat("\n==================================================")
cat("\nDeleting non-symmetrical matrices... ")
non_sym_matrices <- non_sym_matrices[!is.na(non_sym_matrices)]
list_matrices <- list_matrices[-non_sym_matrices]
cat("\nNon-symmetrical matrices deleted!")
# free memory
rm(list=c("list_matrices_mtx", "i"))
invisible(gc()) #called after rm as suggested in R documentation
cat("\n==================================================")
cat("\nMatrices reading finished!")

# Matrices Processing (solving) --------------------------------------------------
cat("\n\nSolving the matrices...")
for(i in seq_along(list_matrices)){
  cat("\n==================================================\n")
  cat("Solving", list_matrices[i], "matrix:")

  # compute B (such that the exact solution of Ax=b is xe=[1,1,..])
  A <- get(list_matrices[i])
  xe <- rep(1,nrow(A))
  b <- crossprod(A, xe)

  # solve Ax=b
  x_benchmark <- mark(x <- SparseM::solve(A, b))

  # time (sec)
  time <- x_benchmark$total_time #time in s
  conv_time <- time/60 #time in min
  cat("\n Execution Time (min): ", conv_time)
  # relative error (norm2)
  relative_error <- norm(x-xe,"2")/norm(xe,"2")
  cat("\n Relative Error (norm2): ", relative_error)
  # RAM used (Bytes)
  conv_ram <- as.numeric(x_benchmark$mem_alloc)
  cat("\n RAM used (Bytes): ", conv_ram)

  # Matrix solving results in dataframe
  matrix_solving_data <- list(df_matrix=list_matrices[i], df_time=conv_time, df_relative_error=relative_error, df_peak_ram=conv_ram)
  matrix_solving_results <- rbind(matrix_solving_results, matrix_solving_data, stringsAsFactors = FALSE)

  # free memory (dealing with big matrix could generate memory errors)
  rm(list=c("A","xe","b","x_benchmark","time", "conv_time", "relative_error","conv_ram",list_matrices[i]))
  invisible(gc()) #called after rm as suggested in R documentation
}
cat("\n==================================================")
cat("\nMatrices solving finished!")

# Matrices POST-Processing (CSV writing) -----------------------------------------
cat("\n\nWriting the results...")
write.csv(matrix_reading_results, file.path(results_dir, "results_matrices_readed.csv"), row.names = FALSE)
write.csv(matrix_solving_results, file.path(results_dir, "results_matrices_solved.csv"), row.names = FALSE)
cat("\nResults written on CSV!")