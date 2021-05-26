# Title     : Direct Methods for Sparse Matrices
# Objective : Solve sparse matrices using an open source programming language (R),
#             a project for the Scientific Computing Methods course of the University of Milan Bicocca
# Created by: Mammarella Davide

# Packages ----------------------------------------------------------------------
cat("\nInitializing the environment...")

# Auto-Mirror + Install only packages not yet installed
suppressMessages({
  chooseCRANmirror(ind=1)
  packages <- c("here","bench","tryCatchLog","Matrix")
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  invisible(lapply(packages, library, character.only = TRUE, quietly=T))
})
rm(list=c("packages", "installed_packages"))

# Initializations ---------------------------------------------------------------
options( warn = -1
         ,tryCatchLog.include.full.call.stack = FALSE
         , tryCatchLog.include.compact.call.stack = FALSE)
set.logging.functions( error.log.func = function(msg) invisible(),
                       warn.log.func  = function(msg) invisible(),
                       info.log.func  = function(msg) invisible())

# Variables initialization
matrices_dir <- "matrices/test/" #matrices directory
list_matrices_mtx <- list.files(path=matrices_dir, pattern=".mtx$", full.names=TRUE) #list of .mtx files in the directory
list_matrices_mtx <- sub(paste0(matrices_dir,"/"),"",
                         list_matrices_mtx[match(seq_along(list_matrices_mtx), rank(file.info(list_matrices_mtx)$size))]) #order by size asc
list_matrices <- sub(".mtx$", "", list_matrices_mtx) #matrices names for loop purpose
results_dir <- "R/results/" #results directory
results_read_csv <- file.path(results_dir, "results_read.csv")
results_solve_csv <- file.path(results_dir, "results_solve.csv")
mtx_read_1row <- data.frame("Matrix Name", "Loading Time (s)", "Size (Bytes)", "Memory Allocation (Bytes)")
mtx_solve_1row <- data.frame("Matrix Name", "Cholesky Solving Time (s)", "x Solving Time (s)", "Total Solving Time (s)", "Relative Error", "Peak RAM (Bytes)")

# CSV Initialization
invisible(file.remove(results_read_csv)) #delete csv files on every start
invisible(file.remove(results_solve_csv)) #delete csv files on every start
write.table(mtx_read_1row, results_read_csv, sep=",", append=TRUE, quote=FALSE, row.names = FALSE, col.names=FALSE)
write.table(mtx_solve_1row, results_solve_csv, sep=",", append=TRUE, quote=FALSE, row.names = FALSE, col.names=FALSE)

# Matrix READING -----------------------------------------------------------------

# Auto-Read from "matrices_dir" + Auto Variable-Name assignation from file (w/o .mtx)
for (i in seq_along(list_matrices)){
  cat("\n==================================================\n")
  cat("Reading", list_matrices[i], "matrix:")
  matrix_name <- sub(".mtx$", "", list_matrices_mtx[i])
  matrix_dir <- paste0(matrices_dir, list_matrices_mtx[i])
  matrix_read_benchmark <- mark(matrix_read <- readMM(matrix_dir), time_unit="s")

  # Symmetric check
  if (isSymmetric(matrix_read)){
    assign(paste0(matrix_name), matrix_read) #assign name to matrix
    matrix_size <- as.numeric(object.size(matrix_read))
    cat("\n Matrix size (Bytes):", matrix_size)
    matrix_allocation <- as.numeric(matrix_read_benchmark$mem_alloc)
    cat("\n Reading memory allocation (Bytes):", matrix_allocation)
    matrix_load_time <- as.numeric(matrix_read_benchmark$total_time) #time in s
    cat("\n Loading time (s):", matrix_load_time)

    # Results in CSV
    row <- data.frame(list_matrices[i], matrix_load_time, matrix_size, matrix_allocation)
    write.table(row, results_read_csv, sep=",", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
    cat("\nResults written in", paste0(results_dir,"results_read.csv"))
  } else {
    cat("\n Matrix non-symmetrical!")
    rm(list=c("matrix_read","matrix_name","matrix_dir"))
    invisible(gc())
    cat("\n Deleted, skipping to next one...")
    next
  }

# Matrices Processing (solving) --------------------------------------------------
  cat("\n==================================================")
  cat("\nSolving", list_matrices[i], "matrix:")

  A <- get(list_matrices[i])

  # detect not-positive definite matrix (cholesky analysis)
  positive_matrix_check <- tryCatchLog({
      R_benchmark <- mark(R <- Cholesky(A))
    }
    , error=function(e){
      cat("\n Matrix not-positive definite!")
      rm(list=c("A","R",list_matrices[i]))
      invisible(gc())
      cat("\n Matrix deleted!")
    }
  )

  # solve positive matrix
  if(!inherits(positive_matrix_check, "error") & exists("A") & exists("R")){
    xe <- rep(1,nrow(A))
    b <- crossprod(A, xe)

    # solve Ax=b
    x_benchmark <- mark(x <- solve(R, b), time_unit="s")

    # time (sec)
    time_chol <- as.numeric(R_benchmark$total_time) #time in s
    conv_time_chol <- time_chol/60 #time in min
    cat("\n Cholesky Resolution Time (min):", conv_time_chol)
    time_x <- as.numeric(x_benchmark$total_time) #time in s
    conv_time_x <- time_x/60 #time in min
    cat("\n x Resolution Time (min):", conv_time_x)
    time_total <- as.numeric(time_x+time_chol)
    conv_time_total <- time_total/60
    cat("\n Total Solving Time (min):", conv_time_total)

    # relative error (norm2)
    relative_error <- norm(x-xe,"2")/norm(xe,"2")
    cat("\n Relative Error (norm2): ", relative_error)

    # RAM used (Bytes)
    # Total amount of memory allocated by R while running the expression.
    # Memory allocated outside the R heap, e.g. by malloc() or new directly is not tracked
    # https://cran.r-project.org/web/packages/bench/bench.pdf
    conv_ram <- as.numeric(x_benchmark$mem_alloc)
    cat("\n RAM used (Bytes): ", conv_ram)

    # Results in CSV
    row1 <- data.frame(list_matrices[i], time_chol, time_x, time_total, relative_error, conv_ram)
    write.table(row1, results_solve_csv, sep=",", append=TRUE, quote=FALSE, row.names = FALSE, col.names=FALSE)
    cat("\nResults written in", paste0(results_dir,"results_solve.csv"))

    # Free MEM
    rm(list=c("x","A","R","xe","b",list_matrices[i]))
    invisible(gc())
  } else {
    # Free MEM
    rm(list=c("A","matrix_read","matrix_name","matrix_dir",list_matrices[i]))
    invisible(gc())
    next
  }
}
cat("\n==================================================")
cat("\nMatrices solving finished!")