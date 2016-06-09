/// \example LOBPCGEpetraExGen.cpp
/// \brief Use LOBPCG with Epetra, for a generalized eigenvalue problem.
///
/// This example computes the eigenvalues of largest magnitude of an
/// generalized eigenvalue problem, using Anasazi's implementation of
/// the LOBPCG method, with Epetra linear algebra.
///
/// The test problem claims to come from research described in the
/// following report: " A comparison of algorithms for modal analysis
/// in the absence of a sparse direct method", P. Arbenz, R. Lehoucq,
/// and U. Hetmaniuk, Sandia National Laboratories, Technical report
/// SAND2003-1028J.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"


using namespace Anasazi;

long int getCurrentMillis()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

int main(int argc, char *argv[]) {
  bool success = false;
  try {
    // Create an Epetra communicator
    //
    Epetra_SerialComm Comm;

    // Create an Anasazi output manager
    //
    BasicOutputManager<double> printer;
    printer.stream(Errors) << Anasazi_Version() << std::endl << std::endl;

    // Get the sorting std::string from the command line
    //
    std::string which = "SM";
    std::string lapFile;
    std::string fvFile;
    int maxIters  = 0;
    double tol = 0;    
    Teuchos::CommandLineProcessor cmdp(false,true); 
    cmdp.setOption("fin",&lapFile, "Filename for Matrix-Market laplacian matrix.");
    cmdp.setOption("fout",&fvFile, "Filename for Fiedler vector.");
    cmdp.setOption("maxIters",&maxIters, "max iters");
    cmdp.setOption("tol",&tol, "tolerance");    
    if(cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
        return -1;
    }
    printer.stream(Errors) << "passed cml ops"  << std::endl;
    printf("maxIters=%d, tol=%3E\n", maxIters, tol);
    
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator OP;
    typedef MultiVecTraits<double, Epetra_MultiVector> MVT;

    // Eigensolver parameters
    int nev = 2;
    int blockSize = 5;
    int errcode = 0;
    
    Epetra_CrsMatrix * Km = NULL;    
    errcode = EpetraExt::MatrixMarketFileToCrsMatrix(lapFile.c_str(),
                                                     Comm,
                                                     Km);
    printer.stream(Errors) << "built Km: " << errcode << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(errcode, std::runtime_error, "error reading file " + lapFile);    
    
    Teuchos::RCP<Epetra_CrsMatrix> K =
        Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(Km), false );
    printer.stream(Errors) << "built K"  << std::endl;    
    
    Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
    ivec->Random();
    printer.stream(Errors) << "passed ivec"  << std::endl;    

    // Create the eigenproblem.
    Teuchos::RCP<BasicEigenproblem<double, MV, OP> > MyProblem =
      Teuchos::rcp( new BasicEigenproblem<double, MV, OP>(K, ivec) );
    printer.stream(Errors) << "passed myprob"  << std::endl;    
    
    // Inform the eigenproblem that the operator A is symmetric
    MyProblem->setHermitian(true);

    // Set the number of eigenvalues requested
    MyProblem->setNEV( nev );

    // Inform the eigenproblem that you are finishing passing it information
    bool boolret = MyProblem->setProblem();
    if (boolret != true) {
      printer.print(Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
      throw -1;
    }
    printer.stream(Errors) << "passed set my prob"  << std::endl;

    // Create parameter list to pass into the solver manager
    //
    Teuchos::ParameterList MyPL;
    MyPL.set( "Which", which );
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Maximum Iterations", maxIters );
    MyPL.set( "Convergence Tolerance", tol );
    MyPL.set( "Full Ortho", true );
    MyPL.set( "Use Locking", true );
    //
    // Create the solver manager
    LOBPCGSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
    printer.stream(Errors) << "passed create solver"  << std::endl;

    // Solve the problem
    //
    long int start = getCurrentMillis();    
    ReturnType returnCode = MySolverMan.solve();
    long int end = getCurrentMillis();
    std::ostringstream os;    
    os << "solver took: " << (end-start) << " ms" << std::endl;
    printer.stream(Errors) << "passed solve"  << std::endl;    

    // Get the eigenvalues and eigenvectors from the eigenproblem
    //
    Eigensolution<double,MV> sol = MyProblem->getSolution();
    std::vector<Value<double> > evals = sol.Evals;
    Teuchos::RCP<MV> evecs = sol.Evecs;
    printer.stream(Errors) << "passed get soln"  << std::endl;     

    // Compute residuals.
    //
    std::vector<double> normR(sol.numVecs);
    double normMat = 0;
    if (sol.numVecs > 0) {
      Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
      Epetra_MultiVector Kvec( K->OperatorDomainMap(), evecs->NumVectors() );
      T.putScalar(0.0);
      for (int i=0; i<sol.numVecs; i++) {
        T(i,i) = evals[i].realpart;
      }
      K->Apply( *evecs, Kvec );
      MVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, Kvec );
      MVT::MvNorm( Kvec, normR );
      normMat = K->NormFrobenius();
    }
    printer.stream(Errors) << "passed calc res"  << std::endl;    

    // Print the results
    //
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os.setf(std::ios::fixed);
    os<<"Solver manager returned " << (returnCode == Converged ? "converged." : "unconverged.") << std::endl;
    os<<std::endl;
    for (int i=0; i<sol.numVecs; i++) {
      printf("ac = %.16f\n", evals[i].realpart);
      printf("relres = %.3E\n", normR[i]/normMat);
    }
    printer.print(Errors,os.str());
    EpetraExt::MultiVectorToMatrixMarketFile(fvFile.c_str(), *evecs);
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
