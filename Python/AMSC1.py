#AMSC1
import scipy as sp
import numpy as np
import pathlib
import time

from scipy.io import mmread
from scipy.sparse import isspmatrix, isspmatrix_csc, isspmatrix_csr
from scipy.sparse.linalg import spsolve
from scipy.linalg import norm
from scipy.sparse import csc_matrix

import psutil #Just for windows


path = str(pathlib.Path(__file__).parent.absolute())

def Matrice_input(path, name_matrix, name_save_file):
    p = psutil.Process()
    start = time.time()
    A = mmread(path + name_matrix)
    end = time.time()
    time_load = end - start
    #print(isspmatrix(A)) #Check is matrix is sparse
    #print(isspmatrix_csc(A))
    #print(isspmatrix_csr(A))
    n_col = A.get_shape()[0]
    x = np.ones(n_col)
    b = A * x
    start = time.time()

    if(isspmatrix_csc(A) or isspmatrix_csr(A) ):
        x_approx = spsolve(A, b)
    else:
        x_approx = spsolve(csc_matrix(A), b)

    end = time.time()
    time_solution = end - start
    err = norm(x - x_approx)/norm(x)

    #For Windows
    #This will return the peak memory usage in bytes.
    ram = psutil.Process().memory_info().peak_wset / 1000000 #MBytes
    #For Linux
    #On Linux this will be measured in KiB
    #ram = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / resource.getpagesize()
     
    #print(ram)
    solution = [[time_load], [time_solution], [err], [ram]]
    a_file = open(path+name_save_file, "w")
    for row in solution:
        np.savetxt(a_file, row)

    a_file.close()

    a_file = open(path+name_save_file[:-4]+"_approx.txt", "w")
    np.savetxt(a_file, x_approx)
    a_file.close()

    return solution

Matrice_input(path, "\\GT01R\\GT01R.mtx", "\\GT01R.txt" ) 
Matrice_input(path, "\\TSC_OPF_1047\\TSC_OPF_1047.mtx", "\\TSC_OPF_1047.txt" )
Matrice_input(path, "\\ns3Da\\ns3Da.mtx", "\\ns3Da.txt" ) 
Matrice_input(path, "\\bundle_adj\\bundle_adj.mtx", "\\bundle_adj.txt" ) 
Matrice_input(path, "\\ifiss_mat\\ifiss_mat.mtx", "\\ifiss_mat.txt" )
Matrice_input(path, "\\G3_circuit\\G3_circuit.mtx", "\\G3_circuit.txt" ) 
Matrice_input(path, "\\nd24k\\nd24k.mtx", "\\nd24k.txt" ) # Va ma lento

Matrice_input(path, "\\Hook_1498\\Hook_1498.mtx", "\\Hook_1498.txt" ) # Can't expand MemType 0: jcol 649435

