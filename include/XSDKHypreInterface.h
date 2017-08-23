// @HEADER
// ***********************************************************************
//
//       xSDKTrilinos: Extreme-scale Software Development Kit Package
//                 Copyright (2016) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alicia Klinvex    (amklinv@sandia.gov)
//                    James Willenbring (jmwille@sandia.gov)
//                    Michael Heroux    (maherou@sandia.gov)         
//
// ***********************************************************************
// @HEADER

// This file was copied into Nalu source tree from
// https://github.com/trilinos/xSDKTrilinos. The contents have been modified to
// enable interfacing with Nalu. The original copyright and license have been
// retained.
//

#ifndef XSDKHYPREINTERFACE_H
#define XSDKHYPREINTERFACE_H

#include "Ifpack2_ConfigDefs.hpp"
//#ifdef HAVE_HYPRE // TODO

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"

#include "Teuchos_RCP.hpp"

#include <cmath>

namespace Ifpack2 {

#ifndef HYPRE_ENUMS
#define HYPRE_ENUMS
namespace Hypre {
//! This enumerated type defines the allowed solvers and preconditioners in Hypre. Some can be used as both solver and preconditioner.
enum Hypre_Solver{
    BoomerAMG,
    ParaSails,
    Euclid,
    AMS,
    Hybrid,
    PCG,
    GMRES,
    FlexGMRES,
    LGMRES,
    BiCGSTAB
};

//! This enumerated type defines the two options for applying inverse, either solve or apply the preconditioner.
enum Hypre_Chooser{
    Solver,
    Prec
};
} // namespace Hypre
#endif //HYPRE_ENUMS

//! This class is used to help with passing parameters in the SetParameter() function. Use this class to call Hypre's internal parameters.
class FunctionParameter{
  public:
    //! Single int constructor.
    FunctionParameter(Hypre::Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int), int param1) :
      chooser_(chooser),
      option_(0),
      int_func_(funct_name),
      int_param1_(param1) {}

    //! Single double constructor.
    FunctionParameter(Hypre::Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, double), double param1):
      chooser_(chooser),
      option_(1),
      double_func_(funct_name),
      double_param1_(param1) {}

    //! Single double, single int constructor.
    FunctionParameter(Hypre::Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, double, int), double param1, int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(funct_name),
      int_param1_(param2),
      double_param1_(param1) {}

    //! Two ints constructor.
    FunctionParameter(Hypre::Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int, int), int param1, int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(funct_name),
      int_param1_(param1),
      int_param2_(param2) {}

    //! Int pointer constructor.
    FunctionParameter(Hypre::Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int*), int *param1):
      chooser_(chooser),
      option_(4),
      int_star_func_(funct_name),
      int_star_param_(param1) {}

    //! Double pointer constructor.
    FunctionParameter(Hypre::Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, double*), double* param1):
      chooser_(chooser),
      option_(5),
      double_star_func_(funct_name),
      double_star_param_(param1) {}

    //! char pointer constructor
    FunctionParameter(Hypre::Hypre_Chooser chooser,
                      int (*funct_name)(HYPRE_Solver, char*),
                      char* param1
    ): chooser_(chooser),
       option_(6),
       char_star_func_(funct_name),
       char_star_param_(param1)
    {}

    //! Only method of this class. Calls the function pointer with the passed in HYPRE_Solver
    int CallFunction(HYPRE_Solver solver, HYPRE_Solver precond){
      if(chooser_ == Hypre::Solver){
        if(option_ == 0){
          return int_func_(solver, int_param1_);
        } else if(option_ == 1){
          return double_func_(solver, double_param1_);
        } else if(option_ == 2){
          return double_int_func_(solver, double_param1_, int_param1_);
        } else if (option_ == 3){
          return int_int_func_(solver, int_param1_, int_param2_);
        } else if (option_ == 4){
          return int_star_func_(solver, int_star_param_);
        } else if (option_ == 5) {
          return double_star_func_(solver, double_star_param_);
        } else if (option_ == 6) {
            return char_star_func_(solver, char_star_param_);
        }
      } else {
        if(option_ == 0){
          return int_func_(precond, int_param1_);
        } else if(option_ == 1){
          return double_func_(precond, double_param1_);
        } else if(option_ == 2){
          return double_int_func_(precond, double_param1_, int_param1_);
        } else if(option_ == 3) {
          return int_int_func_(precond, int_param1_, int_param2_);
        } else if(option_ == 4) {
          return int_star_func_(precond, int_star_param_);
        } else if (option_ == 5) {
            return double_star_func_(precond, double_star_param_);
        } else if (option_ == 6) {
            return char_star_func_(solver, char_star_param_);
        }
      }
      // Should never reach here, but quell warning
      return 0;
    }

  private:
    Hypre::Hypre_Chooser chooser_;
    int option_;
    int (*int_func_)(HYPRE_Solver, int);
    int (*double_func_)(HYPRE_Solver, double);
    int (*double_int_func_)(HYPRE_Solver, double, int);
    int (*int_int_func_)(HYPRE_Solver, int, int);
    int (*int_star_func_)(HYPRE_Solver, int*);
    int (*double_star_func_)(HYPRE_Solver, double*);
    int (*char_star_func_)(HYPRE_Solver, char*);
    int int_param1_;
    int int_param2_;
    double double_param1_;
    int *int_star_param_;
    double *double_star_param_;
    char* char_star_param_;
};

//! Ifpack2_Hypre: A class for constructing and using an ILU factorization of a given Tpetra::RowMatrix, using the Hypre library by Lawrence Livermore National Laboratories.

/*!
Class Ifpack2_Hypre: A class for using methods of Hypre with Tpetra objects.
*/

/*template<class MatrixType, 
         class LocalInverseType = 
         Preconditioner<typename MatrixType::scalar_type,
                        typename MatrixType::local_ordinal_type,
                        typename MatrixType::global_ordinal_type,
                        typename MatrixType::node_type> >
class Ifpack2_Hypre: 
    virtual public Preconditioner<typename MatrixType::scalar_type,
                                  typename MatrixType::local_ordinal_type,
                                  typename MatrixType::global_ordinal_type,
                                  typename MatrixType::node_type> */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class Ifpack2_Hypre: 
    public Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
public:
  // @{ Public Types

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  // @}

  // @{ Constructors and destructors.
  //! Constructor
  Ifpack2_Hypre(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A);

  //! Destructor
  virtual ~Ifpack2_Hypre(){ Destroy();}

  // @}
  // @{ Construction methods

  //! initialize the preconditioner, does not touch matrix values.
  virtual void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool isInitialized() const{ return(isInitialized_);}

  //! Compute ILU factors L and U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the ILU(k) factors.
   */
  virtual void compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool isComputed() const{ return(isComputed_);}


  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes six parameter names: Solver,
     Preconditioner, SolveOrPrecondition, SetPreconditioner, NumFunctions and Functions. These names are
     case sensitive. Solver requires an enumerated parameter of type Hypre_Solver. Preconditioner is similar
     except requires the type be a preconditioner. The options are listed below:
                       Solvers                            Preconditioners
                       BoomerAMG                          BoomerAMG
                       AMS                                ParaSails
                       Hybrid                             AMS
                       PCG (Default)                      Euclid (Default)
                       GMRES
                       FlexGMRES
                       LGMRES
                       BiCGSTAB
     SolveOrPrecondition takes enumerated type Hypre_Chooser, Solver will solve the system, Preconditioner will apply the preconditioner.
     SetPreconditioner takes a boolean, true means the solver will use the preconditioner.
     NumFunctions takes an int that describes how many parameters will be passed into Functions. (This needs to be correct.)
     Functions takes an array of Ref Counted Pointers to an object called FunctionParameter. This class is implemented in Ifpack2_Hypre.h.
     The object takes whether it is Solver or Preconditioner that we are setting a parameter for.
     The function in Hypre that sets the parameter, and the parameters for that function. An example is below:

     RCP<FunctionParameter> functs[2];
     functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000)); // max iterations
     functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, 1e-7)); // conv. tolerance
     list.set("NumFunctions", 2);
     list.set<RCP<FunctionParameter>*>("Functions", functs);
     NOTE: SetParameters() must be called to use ApplyInverse(), the solvers will not be created otherwise. An empty list is acceptable to use defaults.
  */
  void setParameters(const Teuchos::ParameterList& parameterlist);

    //! Set a parameter that takes a single int.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set maximum iterations would be &HYPRE_BoomerAMGSetMaxIter
    \param parameter (In) -The integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int), int parameter);

    //! Set a parameter that takes a single double.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set tolerance would be &HYPRE_BoomerAMGSetTol
    \param parameter (In) -The double parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double), double parameter);

    //! Set a parameter that takes a double then an int.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation weight for a given level would be &HYPRE_BoomerAMGSetLevelRelaxWt
    \param parameter1 (In) -The double parameter being set.
    \param parameter2 (In) - The integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double, int), double parameter1, int parameter2);

    //! Set a parameter that takes two int parameters.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation type for a given level would be &HYPRE_BoomerAMGSetCycleRelaxType
    \param parameter1 (In) -The first integer parameter being set.
    \param parameter2 (In) - The second integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int, int), int parameter1, int parameter2);

    //! Set a parameter that takes a double*.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation weight would be &HYPRE_BoomerAMGSetRelaxWeight
    \param parameter (In) -The double* parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double*), double* parameter);

    //! Set a parameter that takes an int*.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set grid relax type would be &HYPRE_BoomerAMGSetGridRelaxType
    \param parameter (In) -The int* parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int*), int* parameter);

    //! Sets the solver that is used by the Solve() and ApplyInverse() methods. Until this is called, the default solver is PCG.
    /*!
    \param chooser (In) - A Hypre_Chooser enumerated type. If Solver, then we are selecting which solver, if Preconditioner, we are choosing which preconditioner to use.
    \param Solver (In) -A Hypre_Solver enumerated type to select the solver or preconditioner. Options for solver are:
    BoomerAMG, AMS, Hybrid, PCG, GMRES, FlexGMRES, LGMRES, and BiCGSTAB. See Hypre Ref Manual for more info on the solvers.
    Options for Preconditioner are: BoomerAMG, ParaSails, Euclid, and AMS.

    \return Integer error code, set to 0 if successful.
  */

    int SetParameter(Hypre::Hypre_Chooser chooser, Hypre::Hypre_Solver solver);

    //! Sets the solver to use the selected preconditioner.
    /*!
    \param UsePreconditioner (In) -A boolean, true use preconditioner, false do not use the supplied preconditioner with the solver.
    The solver and preconditioner must have been selected and the solver must be one of the following solvers:
      Hybrid, PCG, GMRES, FlexGMRES, LGMRES, BiCGSTAB.

    \return Integer error code, set to 0 if successful.
  */

    int SetParameter(bool UsePreconditioner){ UsePreconditioner_ = UsePreconditioner; return 0;}

    //! Choose to solve the problem or apply the preconditioner.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type, either Solver or Preconditioner.
    The chosen type must have been selected before this method is called.

    \return Integer error code, set to 0 if successful.
  */
    int SetParameter(Hypre::Hypre_Chooser chooser) { SolveOrPrec_ = chooser; return 0;}

  //! Call all the function pointers stored in this object.
    int CallFunctions() const;

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
      does not support transpose use, this method should return a value of -1.

      \param
       UseTranspose_in - (In) If true, multiply by the transpose of operator, otherwise just use operator.

      \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};

  // @}

  // @{ Mathematical functions.

  //! Returns the result of a Tpetra_Operator inverse applied to an Tpetra_MultiVector X in Y.
  /*! In this implementation, we use several existing attributes to determine how virtual
      method apply() should call the concrete method Solve().  We pass in the UpperTriangular(),
      the Tpetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
      if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
      be accounted for when doing a triangular solve.

    \param
           X - (In) A Tpetra_MultiVector of dimension NumVectors to solve for.
    \param Out
           Y - (Out) A Tpetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */

  virtual void apply(const Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >& X,
             Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Computes the estimated condition number and returns the value.
  double Condest(const Ifpack2::CondestType CT = Ifpack2::Cheap,
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
                 Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in = 0);

  //! Returns the computed estimated condition number, or -1.0 if not computed.
  double Condest() const{ return(Condest_);}

  // @}
  // @{ Query methods

  //! Returns a character string describing the operator
  const char* Label() const {return(Label_);}

  //! Sets label for \c this object.
  int SetLabel(const char* Label_in)
  {
    strcpy(Label_,Label_in);
    return(0);
  }

  //! Returns a reference to the map that should be used for domain.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const{ return A_->getDomainMap();}

  //! Returns a reference to the map that should be used for range.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const{ return A_->getRangeMap();}

  //! Returns 0.0 because this class cannot compute Inf-norm.
  double NormInf() const {return(0.0);};

  //! Returns false because this class cannot compute an Inf-norm.
  bool HasNormInf() const {return(false);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns the Tpetra_BlockMap object associated with the range of this matrix operator.
  Teuchos::RCP<const Teuchos::Comm<int> > Comm() const{return(A_->getComm());};

  //! Returns a reference to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getMatrix() const{ return(A_);}

  //! Returns the Hypre matrix that was created upon construction.
  const HYPRE_IJMatrix& HypreMatrix()
  {
    if(isInitialized() == false)
      initialize();
    return(HypreA_);
  }

  //! Prints on stream basic information about \c this object.
  virtual std::ostream& Print(std::ostream& os) const;

  //! Returns the number of calls to initialize().
  virtual int getNumInitialize() const{ return(numInitialize_);}

  //! Returns the number of calls to compute().
  virtual int getNumCompute() const{ return(numCompute_);}

  //! Returns the number of calls to ApplyInverse().
  virtual int getNumApply() const{ return(numApply_);}

  //! Returns the time spent in initialize().
  virtual double getInitializeTime() const{ return(initializeTime_);}

  //! Returns the time spent in compute().
  virtual double getComputeTime() const{ return(computeTime_);}

  //! Returns the time spent in ApplyInverse().
  virtual double getApplyTime() const{ return(applyTime_);}

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const{ return(0.0);}

  //! Returns the number of flops in the compute phase.
  virtual double ComputeFlops() const{ return(ComputeFlops_);}

  //! Returns the number of flops in the apply inverse phase.
  virtual double ApplyInverseFlops() const{ return(ApplyInverseFlops_);}

    //! Destroys all internal data
    virtual void Destroy();


    //! Number of linear iterations
    mutable int hypreNumLinIters_;

    //! Final residuals
    mutable double hypreFinalResidual_;
protected:

  // @}
  // @{ Private methods

  //! Copy constructor (should never be used)
  Ifpack2_Hypre(const Ifpack2_Hypre& RHS) {}

  //! operator= (should never be used)
  Ifpack2_Hypre& operator=(const Ifpack2_Hypre& RHS){ return(*this);}

  //! Returns the MPI communicator used in the Tpetra matrix
  MPI_Comm GetMpiComm() const
    { return *(Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A_->getComm()))->getRawMpiComm(); }

  //! Returns the result of a Ifpack_ILU forward/back solve on a Tpetra_MultiVector X in Y.
  /*!
    \param In
    Trans -If true, solve transpose problem.
    \param
    X - (In) A Tpetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y - (Out) A Tpetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >& X, Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >& Y) const;


  //! Returns the number of global matrix rows.
  int getGlobalNumRows() const {return(A_->getGlobalNumRows());};

  //! Returns the number of global matrix columns.
  int getGlobalNumCols() const {return(A_->getGlobalNumCols());};

  //! Returns the number of local matrix rows.
  int getNodeNumRows() const {return(A_->getNodeNumRows());};

  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(A_->NumMyCols());};

  //! Sets the solver type to be the passed in solver type.
  int SetSolverType(Hypre::Hypre_Solver solver);

  //! Sets the preconditioner type to be the passed in type.
  int SetPrecondType(Hypre::Hypre_Solver precond);

  //! Create the solver.
  int CreateSolver();

  //! Create the Preconditioner.
  int CreatePrecond();

  //! Add a function to be called in compute()
  int AddFunToList(Teuchos::RCP<FunctionParameter> NewFun);

  //! Create a BoomerAMG solver.
  int Hypre_BoomerAMGCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_BoomerAMGCreate(solver);}

  //! Create a ParaSails solver.
  int Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParaSailsCreate(comm, solver);}

  //! Create a Euclid solver.
  int Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_EuclidCreate(comm, solver);}

  //! Create an AMS solver.
  int Hypre_AMSCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_AMSCreate(solver);}

  //! Create a Hybrid solver.
  int Hypre_ParCSRHybridCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRHybridCreate(solver);}

  //! Create a PCG solver.
  int Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRPCGCreate(comm, solver);}

  //! Create a GMRES solver.
  int Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRGMRESCreate(comm, solver);}

  //! Create a FlexGMRES solver.
  int Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRFlexGMRESCreate(comm, solver);}

  //! Create a LGMRES solver.
  int Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRLGMRESCreate(comm, solver);}

  //! Create a BiCGSTAB solver.
  int Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRBiCGSTABCreate(comm, solver);}

    void processLocalMatrixInfo();

    void initHypreDataStructures();

  // @}
  // @{ Internal data

  //! Pointer to the Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> to factorize
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! This objects copy of the ParameterList
  Teuchos::ParameterList List_;
  //! Needed to support Tpetra_Operator abstract class
  bool UseTranspose_;
  //! A condition estimate for the preconditioner, -1 until compute()
  double Condest_;
  //! If \c true, the preconditioner has been successfully initialized.
  bool isInitialized_;
  //! If \c true, the preconditioner has been successfully computed.
  bool isComputed_;
  //! Label of \c this object.
  char Label_[160];
  //! Contains the number of successful calls to initialize().
  int numInitialize_;
  //! Contains the number of successful call to compute().
  int numCompute_;
  //! Contains the number of successful call to apply().
  mutable int numApply_;
  //! Contains the time for all successful calls to initialize().
  double initializeTime_;
  //! Contains the time for all successful calls to compute().
  double computeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double applyTime_;
  //! Contains the number of flops for compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;

  //! The Hypre matrix created in initialize()
  mutable HYPRE_IJMatrix HypreA_;
  //! Pointer to the CSR (same matrix)
  mutable HYPRE_ParCSRMatrix ParMatrix_;
  //! The Hypre Vector for input
  mutable HYPRE_IJVector XHypre_;
  //! The Hypre Vector for output
  mutable HYPRE_IJVector YHypre_;
  mutable HYPRE_ParVector ParX_;
  mutable HYPRE_ParVector ParY_;
  mutable hypre_ParVector *XVec_;
  mutable hypre_ParVector *YVec_;
  mutable hypre_Vector *XLocal_;
  mutable hypre_Vector *YLocal_;
  //! The Hypre Solver if doing a solve
  mutable HYPRE_Solver Solver_;
  //! The Hypre Solver if applying preconditioner
  mutable HYPRE_Solver Preconditioner_;
  //  The following are pointers to functions to use the solver and preconditioner.
  int (Ifpack2_Hypre::*SolverCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  int (*SolverDestroyPtr_)(HYPRE_Solver);
  int (*SolverSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*SolverSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*SolverPrecondPtr_)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  int (Ifpack2_Hypre::*PrecondCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  int (*PrecondDestroyPtr_)(HYPRE_Solver);
  int (*PrecondSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*PrecondSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);

    int (*SolverNumItersPtr_)(HYPRE_Solver, int*);
    int (*SolverFinalResidualNormPtr_)(HYPRE_Solver, double*);

  bool *IsSolverSetup_;
  bool *IsPrecondSetup_;
  //! Is the system to be solved or apply preconditioner
  Hypre::Hypre_Chooser SolveOrPrec_;
  //! Counter of the number of parameters set
  int NumFunsToCall_;
  //! Which solver was chosen
  Hypre::Hypre_Solver SolverType_;
  //! Which preconditioner was chosen
  Hypre::Hypre_Solver PrecondType_;
  //! Should the preconditioner be used in the solver
  bool UsePreconditioner_;
  //! This contains a list of function pointers that will be called in compute
  std::vector<Teuchos::RCP<FunctionParameter> > FunsToCall_;

    int hypreILower_;
    int hypreIUpper_;

};

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2_Hypre<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initHypreDataStructures()
{
    MPI_Comm comm = GetMpiComm();
    // Next create vectors that will be used when ApplyInverse() is called
    HYPRE_IJVectorCreate(comm, hypreILower_, hypreIUpper_, &XHypre_);
    HYPRE_IJVectorSetObjectType(XHypre_, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(XHypre_);
    HYPRE_IJVectorAssemble(XHypre_);
    HYPRE_IJVectorGetObject(XHypre_, (void**) &ParX_);

    HYPRE_IJVectorCreate(comm, hypreILower_, hypreIUpper_, &YHypre_);
    HYPRE_IJVectorSetObjectType(YHypre_, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(YHypre_);
    HYPRE_IJVectorAssemble(YHypre_);
    HYPRE_IJVectorGetObject(YHypre_, (void**) &ParY_);

    XVec_ = (hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) XHypre_));
    XLocal_ = hypre_ParVectorLocalVector(XVec_);

    YVec_ = (hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) YHypre_));
    YLocal_ = hypre_ParVectorLocalVector(YVec_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::processLocalMatrixInfo()
{
    MPI_Comm comm = GetMpiComm();
    int nprocs, iproc;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &iproc);
    int nLocalRows = A_->getNodeNumRows();

    std::vector<int> indexMap(nprocs);
    std::vector<int> startIdx(nprocs+1);
    MPI_Allgather(&nLocalRows, 1, MPI_INT, indexMap.data(), 1, MPI_INT, comm);

    startIdx[0] = 0;
    for (int i=1; i<=nprocs; i++)
        startIdx[i] = startIdx[i-1] + indexMap[i-1];

    hypreILower_ = startIdx[iproc] ;
    hypreIUpper_ = hypreILower_ + nLocalRows - 1;
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Ifpack2_Hypre(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A):
  A_(A),
  UseTranspose_(false),
  isInitialized_(false),
  isComputed_(false),
  Label_(),
  numInitialize_(0),
  numCompute_(0),
  numApply_(0),
  initializeTime_(0.0),
  computeTime_(0.0),
  applyTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  SolveOrPrec_(Hypre::Solver),
  NumFunsToCall_(0),
  SolverType_(Hypre::PCG),
  PrecondType_(Hypre::Euclid),
  UsePreconditioner_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!A_->isFillComplete(),std::invalid_argument,
      "Ifpack2::Hypre: Please call fillComplete and try again.");

  IsSolverSetup_ = new bool[1];
  IsPrecondSetup_ = new bool[1];
  IsSolverSetup_[0] = false;
  IsPrecondSetup_[0] = false;

  TEUCHOS_TEST_FOR_EXCEPTION(!A_->getDomainMap()->isSameAs(*A_->getRangeMap()), std::runtime_error,
      Teuchos::typeName (*this) << "::apply(): A's domain and range map must be the same for hypre.");

} //Constructor

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Destroy(){
  if(isInitialized()){
    HYPRE_IJMatrixDestroy(HypreA_);
  }
  HYPRE_IJVectorDestroy(XHypre_);
  HYPRE_IJVectorDestroy(YHypre_);
  if(IsSolverSetup_[0]){
    SolverDestroyPtr_(Solver_);
  }
  if(IsPrecondSetup_[0]){
    PrecondDestroyPtr_(Preconditioner_);
  }
  delete[] IsSolverSetup_;
  delete[] IsPrecondSetup_;
} //Destroy()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initialize(){
  using Teuchos::RCP;

  // Create a timer
  const std::string timerName ("Ifpack2::Hypre::initialize");
  RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  { // Start timer here
    Teuchos::TimeMonitor timeMon (*timer);

    MPI_Comm comm = GetMpiComm();

    int nprocs, iproc;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &iproc);

    HYPRE_IJMatrixCreate(comm, hypreILower_, hypreIUpper_, hypreILower_, hypreIUpper_, &HypreA_);
    HYPRE_IJMatrixSetObjectType(HypreA_, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(HypreA_);

    for(size_t i = 0; i < A_->getNodeNumRows(); i++){
      int numElements = A_->getNumEntriesInLocalRow(i);
      Teuchos::Array<LocalOrdinal> indices(numElements);
      Teuchos::Array<Scalar> values(numElements);
      size_t numEntries;
      A_->getLocalRowCopy(i, indices(), values(), numEntries);
      for(size_t j = 0; j < numEntries; j++){
          int gidx = A_->getColMap()->getGlobalElement(indices[j]) -1;
           indices[j] = gidx;
      }
      int GlobalRow[1];
      GlobalRow[0] = A_->getRowMap()->getGlobalElement(i) - 1;
      HYPRE_IJMatrixAddToValues(HypreA_, 1, &numElements, GlobalRow, &indices[0], &values[0]);
    }
    HYPRE_IJMatrixAssemble(HypreA_);
    HYPRE_IJMatrixGetObject(HypreA_, (void**)&ParMatrix_);

  } // Stop timer here

  isInitialized_=true;
  numInitialize_++;
  initializeTime_ = timer->totalElapsedTime();
} //initialize()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameters(const Teuchos::ParameterList& list){
  using Teuchos::RCP;

  List_ = list;
  Hypre::Hypre_Solver solType = List_.get("Solver", Hypre::PCG);
  SolverType_ = solType;
  Hypre::Hypre_Solver precType = List_.get("Preconditioner", Hypre::Euclid);
  PrecondType_ = precType;
  Hypre::Hypre_Chooser chooser = List_.get("SolveOrPrecondition", Hypre::Solver);
  SolveOrPrec_ = chooser;
  bool SetPrecond = List_.get("SetPreconditioner", false);
  SetParameter(SetPrecond);
  int NumFunctions = List_.get("NumFunctions", 0);
  FunsToCall_.clear();
  NumFunsToCall_ = 0;
  if(NumFunctions > 0){
    RCP<FunctionParameter>* params = List_.get<RCP<FunctionParameter>*>("Functions");
    for(int i = 0; i < NumFunctions; i++){
      AddFunToList(params[i]);
    }
  }
} //SetParameters()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::AddFunToList(Teuchos::RCP<FunctionParameter> NewFun){
  NumFunsToCall_ = NumFunsToCall_+1;
  FunsToCall_.resize(NumFunsToCall_);
  FunsToCall_[NumFunsToCall_-1] = NewFun;
  return 0;
} //AddFunToList()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int), int parameter){
  Teuchos::RCP<FunctionParameter> temp = Teuchos::rcp(new FunctionParameter(chooser, pt2Func, parameter));
  AddFunToList(temp);
  return 0;
} //SetParameter() - int function pointer

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double), double parameter){
  Teuchos::RCP<FunctionParameter> temp = Teuchos::rcp(new FunctionParameter(chooser, pt2Func, parameter));
  AddFunToList(temp);
  return 0;
} //SetParameter() - double function pointer

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double, int), double parameter1, int parameter2){
  Teuchos::RCP<FunctionParameter> temp = Teuchos::rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  AddFunToList(temp);
  return 0;
} //SetParameter() - double,int function pointer

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int, int), int parameter1, int parameter2){
  Teuchos::RCP<FunctionParameter> temp = Teuchos::rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  AddFunToList(temp);
  return 0;
} //SetParameter() int,int function pointer

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double*), double* parameter){
  Teuchos::RCP<FunctionParameter> temp = Teuchos::rcp(new FunctionParameter(chooser, pt2Func, parameter));
  AddFunToList(temp);
  return 0;
} //SetParameter() - double* function pointer

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameter(Hypre::Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int*), int* parameter){
  Teuchos::RCP<FunctionParameter> temp = Teuchos::rcp(new FunctionParameter(chooser, pt2Func, parameter));
  AddFunToList(temp);
  return 0;
} //SetParameter() - int* function pointer

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameter(Hypre::Hypre_Chooser chooser, Hypre::Hypre_Solver solver){
  if(chooser == Hypre::Solver){
    SolverType_ = solver;
  } else {
    PrecondType_ = solver;
  }
  return 0;
} //SetParameter() - set type of solver

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::compute(){
  // Create a timer
  const std::string timerName ("Ifpack2::Hypre::compute");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  if(isInitialized() == false){
    initialize();
  }

  { // Start timer here
    Teuchos::TimeMonitor timeMon (*timer);

    SetSolverType(SolverType_);
    SetPrecondType(PrecondType_);
    CallFunctions();
    if(UsePreconditioner_){
      if(SolverPrecondPtr_ != NULL){
        SolverPrecondPtr_(Solver_, PrecondSolvePtr_, PrecondSetupPtr_, Preconditioner_);
      }
    }
    if(SolveOrPrec_ == Hypre::Solver){
      SolverSetupPtr_(Solver_, ParMatrix_, ParX_, ParY_);
      IsSolverSetup_[0] = true;
    } else {
      PrecondSetupPtr_(Preconditioner_, ParMatrix_, ParX_, ParY_);
      IsPrecondSetup_[0] = true;
    }
  } // Stop timer here

  isComputed_ = true;
  numCompute_++;
  computeTime_ = timer->totalElapsedTime();
} //compute()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CallFunctions() const{
  for(int i = 0; i < NumFunsToCall_; i++){
    FunsToCall_[i]->CallFunction(Solver_, Preconditioner_);
  }
  return 0;
} //CallFunctions()

//==============================================================================
// TODO: Unlike the Tpetra interface, I assume that X != Y.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(const Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >& X, 
                         Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >& Y,
                         Teuchos::ETransp mode, // TODO: These parameters currently do nothing
                         Scalar alpha,
                         Scalar beta) const
{
  // Create a timer
  const std::string timerName ("Ifpack2::Hypre::apply");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
         Teuchos::typeName (*this) << "::apply(): Preconditioner has not been computed.");

  TEUCHOS_TEST_FOR_EXCEPTION(!A_->getDomainMap()->isSameAs(*X.getMap()), std::runtime_error,
      Teuchos::typeName (*this) << "::apply(): X's map must match A's domain map.");

  TEUCHOS_TEST_FOR_EXCEPTION(!A_->getRangeMap()->isSameAs(*Y.getMap()), std::runtime_error,
      Teuchos::typeName (*this) << "::apply(): Y's map must match A's range map.");

  { // Start timer here
    Teuchos::TimeMonitor timeMon (*timer);

    size_t NumVectors = X.getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(NumVectors != Y.getNumVectors(), std::runtime_error,
         Teuchos::typeName (*this) << "::apply(): X and Y must have the same number of vectors.");

    for (size_t VecNum=0; VecNum < NumVectors; VecNum++) {
        auto rhsVec = X.getVector(VecNum);
        auto rhsView = rhsVec->get1dView();
        auto rhsMap = rhsVec->getMap();
        auto limin = rhsMap->getMinLocalIndex();
        auto limax = rhsMap->getMaxLocalIndex();

        for (int i=limin; i <= limax; i++) {
            int ig = rhsMap->getGlobalElement(i) - 1;
            double value = rhsView[i];
            HYPRE_IJVectorAddToValues(XHypre_, 1, &ig, &value);
        }
        HYPRE_IJVectorAssemble(XHypre_);
        HYPRE_IJVectorGetObject(XHypre_, (void**) &ParX_);

        HYPRE_ParVectorSetConstantValues(ParY_, 0.0);
        if(SolveOrPrec_ == Hypre::Solver){
            // Use the solver methods
            SolverSolvePtr_(Solver_, ParMatrix_, ParX_, ParY_);
        } else {
            // Apply the preconditioner
            PrecondSolvePtr_(Preconditioner_, ParMatrix_, ParX_, ParY_);
        }

        auto slnData = Y.getDataNonConst(VecNum);
        double outval;
        for (int i=limin; i <= limax; i++) {
            int ig = rhsMap->getGlobalElement(i) - 1;
            HYPRE_IJVectorGetValues(YHypre_, 1, &ig, &outval);
            slnData[i] = outval;
        }
        SolverNumItersPtr_(Solver_, &hypreNumLinIters_);
        SolverFinalResidualNormPtr_(Solver_, &hypreFinalResidual_);
    }
  } // Stop timer here

  numApply_++;
  applyTime_ = timer->totalElapsedTime();
} //ApplyInverse()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::ostream& Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Print(std::ostream& os) const{
  using std::endl;

  if (!Comm()->getRank()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack2_Hypre: " << Label () << endl << endl;
    os << "Using " << Comm()->getSize() << " processors." << endl;
    os << "Global number of rows            = " << A_->getGlobalNumRows() << endl;
    os << "Global number of nonzeros        = " << A_->getGlobalNumEntries() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "initialize()    "   << std::setw(5) << numInitialize_
       << "  " << std::setw(15) << initializeTime_
       << "              0.0              0.0" << endl;
    os << "compute()       "   << std::setw(5) << numCompute_
       << "  " << std::setw(15) << computeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (computeTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / computeTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << numApply_
       << "  " << std::setw(15) << applyTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (applyTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / applyTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }
  return os;
} //Print()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Condest(const Ifpack2::CondestType CT,
                             const int MaxIters,
                             const double Tol,
                             Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in){
  if (!isComputed()) // cannot compute right now
    return(-1.0);
  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);
  return(Condest_);
} //Condest()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetSolverType(Hypre::Hypre_Solver solver){
  switch(solver) {
    case Hypre::BoomerAMG:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_BoomerAMGCreate;
      SolverDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      SolverSetupPtr_ = &HYPRE_BoomerAMGSetup;
      SolverPrecondPtr_ = NULL;
      SolverSolvePtr_ = &HYPRE_BoomerAMGSolve;
      SolverNumItersPtr_ = &HYPRE_BoomerAMGGetNumIterations;
      SolverFinalResidualNormPtr_ = &HYPRE_BoomerAMGGetFinalRelativeResidualNorm;
      break;
    case Hypre::AMS:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_AMSCreate;
      SolverDestroyPtr_ = &HYPRE_AMSDestroy;
      SolverSetupPtr_ = &HYPRE_AMSSetup;
      SolverSolvePtr_ = &HYPRE_AMSSolve;
      SolverPrecondPtr_ = NULL;
      SolverNumItersPtr_ = &HYPRE_AMSGetNumIterations;
      SolverFinalResidualNormPtr_ = &HYPRE_AMSGetFinalRelativeResidualNorm;
      break;
    case Hypre::Hybrid:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_ParCSRHybridCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRHybridDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRHybridSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRHybridSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRHybridSetPrecond;
      break;
    case Hypre::PCG:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_ParCSRPCGCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRPCGDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRPCGSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRPCGSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRPCGSetPrecond;
      SolverNumItersPtr_ = &HYPRE_PCGGetNumIterations;
      SolverFinalResidualNormPtr_ = &HYPRE_PCGGetFinalRelativeResidualNorm;
      break;
    case Hypre::GMRES:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_ParCSRGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRGMRESSetPrecond;
      SolverNumItersPtr_ = &HYPRE_GMRESGetNumIterations;
      SolverFinalResidualNormPtr_ = &HYPRE_GMRESGetFinalRelativeResidualNorm;
      break;
    case Hypre::FlexGMRES:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_ParCSRFlexGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRFlexGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRFlexGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRFlexGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRFlexGMRESSetPrecond;
      SolverNumItersPtr_ = &HYPRE_FlexGMRESGetNumIterations;
      SolverFinalResidualNormPtr_ = &HYPRE_FlexGMRESGetFinalRelativeResidualNorm;
      break;
    case Hypre::LGMRES:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_ParCSRLGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRLGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRLGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRLGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRLGMRESSetPrecond;
      break;
    case Hypre::BiCGSTAB:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack2_Hypre::Hypre_ParCSRBiCGSTABCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRBiCGSTABDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRBiCGSTABSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRBiCGSTABSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRBiCGSTABSetPrecond;
      break;
    default:
      return -1;
    }
  CreateSolver();
  return 0;
} //SetSolverType()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetPrecondType(Hypre::Hypre_Solver Precond){
  switch(Precond) {
    case Hypre::BoomerAMG:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack2_Hypre::Hypre_BoomerAMGCreate;
      PrecondDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      PrecondSetupPtr_ = &HYPRE_BoomerAMGSetup;
      PrecondSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case Hypre::ParaSails:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack2_Hypre::Hypre_ParaSailsCreate;
      PrecondDestroyPtr_ = &HYPRE_ParaSailsDestroy;
      PrecondSetupPtr_ = &HYPRE_ParaSailsSetup;
      PrecondSolvePtr_ = &HYPRE_ParaSailsSolve;
      break;
    case Hypre::Euclid:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack2_Hypre::Hypre_EuclidCreate;
      PrecondDestroyPtr_ = &HYPRE_EuclidDestroy;
      PrecondSetupPtr_ = &HYPRE_EuclidSetup;
      PrecondSolvePtr_ = &HYPRE_EuclidSolve;
      break;
    case Hypre::AMS:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack2_Hypre::Hypre_AMSCreate;
      PrecondDestroyPtr_ = &HYPRE_AMSDestroy;
      PrecondSetupPtr_ = &HYPRE_AMSSetup;
      PrecondSolvePtr_ = &HYPRE_AMSSolve;
      break;
    default:
      return -1;
    }
  CreatePrecond();
  return 0;

} //SetPrecondType()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CreateSolver(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*SolverCreatePtr_)(comm, &Solver_);
} //CreateSolver()

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Ifpack2_Hypre<Scalar,LocalOrdinal,GlobalOrdinal,Node>::CreatePrecond(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*PrecondCreatePtr_)(comm, &Preconditioner_);
} //CreatePrecond()

} // namespace Ifpack2

//#endif // HAVE_HYPRE
#endif /* XSDKHYPREINTERFACE_H */
