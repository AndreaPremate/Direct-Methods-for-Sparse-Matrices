#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
//#include <Eigen/SparseLU>
#include <chrono>
//#include <Eigen/SparseCholesky>
#include <Eigen/Core>
#include <fstream>
#include <string>
using Eigen::MatrixXd;

//windows
#include <windows.h>
#include <stdio.h>
#include <psapi.h>

int PrintMemoryInfo( DWORD processID ){

    HANDLE hProcess;
    PROCESS_MEMORY_COUNTERS pmc;

    // Print the process identifier.
    //printf( "\nProcess ID: %u\n", processID );
    // Print information about the memory usage of the process.

    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
                                    PROCESS_VM_READ,
                                    FALSE, processID );
    unsigned int a = 1000000;
    if (NULL == hProcess)
        return -1;

    if ( GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc)) ){

        /*lasciare commentato
        printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
        
        printf( "\tPeakWorkingSetSize: %u \n", // 0x%08X esadecimale
                  (pmc.PeakWorkingSetSize/a));
        
        printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
        printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n", 
                  pmc.QuotaPeakPagedPoolUsage );
        printf( "\tQuotaPagedPoolUsage: 0x%08X\n", 
                  pmc.QuotaPagedPoolUsage );
        printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n", 
                  pmc.QuotaPeakNonPagedPoolUsage );
        printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n", 
                  pmc.QuotaNonPagedPoolUsage );
        printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage ); 
        printf( "\tPeakPagefileUsage: 0x%08X\n", 
                  pmc.PeakPagefileUsage );
         */      	   
    }

    CloseHandle( hProcess );
    return pmc.PeakWorkingSetSize/a;
}
//fine_windows


/*//  LINUX MEMORY
#include <unistd.h>
#include <sys/resource.h>
//fine linux*/


void solve_CHOLESKY(const char *__s){
    //MATRIX AND XE LOADING

    std::ofstream filewrite;
    std::string res = " Results_Cholesky.txt";
    filewrite.open(__s + res);

    typedef Eigen::SparseMatrix<double> SMatrixXd;
    SMatrixXd A;
    std::chrono::steady_clock::time_point begin_load = std::chrono::steady_clock::now();             //start load
    Eigen::loadMarket(A, __s);
    Eigen::SparseMatrix<double> Resulting_Matrix = A.selfadjointView<Eigen::Lower>();
    //Eigen::SparseMatrix<double> result2 = A;
    std::chrono::steady_clock::time_point end_load = std::chrono::steady_clock::now();               //end load
    filewrite << "Time spent loading the matrix = " << std::chrono::duration_cast<std::chrono::seconds>(end_load - begin_load).count() << "[s]" << std::endl;
    Resulting_Matrix.makeCompressed();
   
    Eigen::VectorXd xe = Eigen::VectorXd::Ones(Resulting_Matrix.cols());
    Eigen::VectorXd b = Resulting_Matrix*xe;

    //CHOLESKY DECOMPOSITION
    std::chrono::steady_clock::time_point begin_LDLT = std::chrono::steady_clock::now();             //start LDLT

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>   solverLDLT;
    solverLDLT.analyzePattern(Resulting_Matrix);
    solverLDLT.factorize(Resulting_Matrix);
    Eigen::VectorXd xLDLT = solverLDLT.solve(b);
    
    /*if(solverLDLT.info() == Eigen::NumericalIssue)
        std::cout << "LDLT SUCCESS" << std::endl;
    */

    std::chrono::steady_clock::time_point end_LDLT = std::chrono::steady_clock::now();               //end LU LDLT

    filewrite << "Time spent using LDLT decomposition = " << std::chrono::duration_cast<std::chrono::seconds>(end_LDLT - begin_LDLT).count() << "[s]" << std::endl;

    double relative_errorLDLT = (xLDLT - xe).norm() / xe.norm();
    filewrite << "The relative error using LDLT decomposition is: " << relative_errorLDLT << std::endl;

    //windows
    filewrite << "PeakWorkingSetSize(windows): " << PrintMemoryInfo(GetCurrentProcessId()) << " MB"<< std::endl;
	 //fine windows 
	
    /*//linux
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    filewrite << "PeakWorkingSetSize(linux): " << (size_t)(rusage.ru_maxrss * 1024L)/1000000 << " MB"<< std::endl;
	//fine linux*/
	
    filewrite.close();

    std::ofstream writeSolution;
    std::string sol = " solution.txt";
    writeSolution.open(__s + sol);
    writeSolution << xLDLT - xe << std::endl;
    writeSolution.close();
}


void solve_LU(const char *__s, std::string symmetric){

    //MATRIX AND XE LOADING

    std::ofstream filewrite;
    std::string res = " Results_LU.txt";
    filewrite.open(__s + res);

    typedef Eigen::SparseMatrix<double> SMatrixXd;
    SMatrixXd A;
    std::chrono::steady_clock::time_point begin_load = std::chrono::steady_clock::now();             //start load
    Eigen::loadMarket(A, __s);
    Eigen::SparseMatrix<double> Resulting_Matrix;
    if(symmetric.compare("non_symmetric") == 0)
        Resulting_Matrix = A;      //PER MATRICI NON SIMMETRICHE
    else if(symmetric.compare("symmetric") == 0)
        Resulting_Matrix = A.selfadjointView<Eigen::Lower>();
    else {
        std::cout << "Second paramater error"<< std::endl;
        return ; 
    }

    std::chrono::steady_clock::time_point end_load = std::chrono::steady_clock::now();               //end load
    filewrite << "Time spent loading the matrix = " << std::chrono::duration_cast<std::chrono::seconds>(end_load - begin_load).count() << "[s]" << std::endl;
    Resulting_Matrix.makeCompressed();

    Eigen::VectorXd xe = Eigen::VectorXd::Ones(Resulting_Matrix.cols());
    Eigen::VectorXd b = Resulting_Matrix*xe;

    //LU DECOMPOSITION

    std::chrono::steady_clock::time_point begin_LU = std::chrono::steady_clock::now();             //start LU

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solverLU; //, Eigen::COLAMDOrdering<int> 
    solverLU.analyzePattern(Resulting_Matrix);
    solverLU.factorize(Resulting_Matrix);
    if (solverLU.info() == Eigen::NumericalIssue)
        std::cout<< solverLU.info() << Eigen::NumericalIssue <<  std::endl;
    Eigen::VectorXd xLU = solverLU.solve(b);

    std::chrono::steady_clock::time_point end_LU = std::chrono::steady_clock::now();               //end LU

    filewrite << "Time spent using LU decomposition = " << std::chrono::duration_cast<std::chrono::seconds>(end_LU - begin_LU).count() << "[s]" << std::endl;

    double relative_errorLU = (xLU - xe).norm() / xe.norm();
    filewrite << "The relative error using LU decomposition is: " << relative_errorLU << std::endl;

     //windows
    filewrite << "PeakWorkingSetSize(windows): " << PrintMemoryInfo(GetCurrentProcessId()) << " MB"<< std::endl;
	 //fine windows 
	 
    /*//linux
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    filewrite << "PeakWorkingSetSize(linux): " << (size_t)(rusage.ru_maxrss * 1024L)/1000000 << " MB"<< std::endl;
	//fine linux*/
	
    filewrite.close();

    std::ofstream writeSolution;
    std::string sol = " solution.txt";
    writeSolution.open(__s + sol);
    writeSolution << xLU - xe << std::endl;
    writeSolution.close();

}


int main(){
    
    //eseguire una matrice alla volta per stabilire il picco di ram esatto(sennò torna solo quello maggiore)

    //solve_LU("GT01R.mtx", "non_symmetric");
    //solve_LU("TSC_OPF_1047.mtx", "symmetric");
    //solve_LU("ns3Da.mtx", "non_symmetric");
    //solve_CHOLESKY("bundle_adj.mtx");
    //solve_LU("ifiss_mat.mtx", "non_symmetric");
    solve_CHOLESKY("G3_circuit.mtx");
    //solve_CHOLESKY("nd24k.mtx");
    //solve_CHOLESKY("Hook_1498.mtx");
}