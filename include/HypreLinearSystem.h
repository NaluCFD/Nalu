/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HYPRELINEARSYSTEM_H
#define HYPRELINEARSYSTEM_H

#include "LinearSystem.h"
#include "overset/OversetManager.h"
#include "overset/OversetInfo.h"

#include "stk_mesh/base/BulkData.hpp"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"
#include "krylov.h"
#include "HYPRE.h"

#include <unordered_set>

namespace sierra {
namespace nalu {

/** Nalu interface to populate a Hypre Linear System
 */
class HypreLinearSystem : public LinearSystem
{
public:
  HypreLinearSystem(
    Realm& realm,
    const unsigned numDof,
    EquationSystem *eqSys,
    LinearSolver *linearSolver);

  virtual ~HypreLinearSystem();

  // Graph/Matrix Construction
  virtual void buildNodeGraph(const stk::mesh::PartVector&);// for nodal assembly (e.g., lumped mass and source)
  virtual void buildFaceToNodeGraph(const stk::mesh::PartVector&);// face->node assembly
  virtual void buildEdgeToNodeGraph(const stk::mesh::PartVector&);// edge->node assembly
  virtual void buildElemToNodeGraph(const stk::mesh::PartVector&);// elem->node assembly
  virtual void buildReducedElemToNodeGraph(const stk::mesh::PartVector&);// elem (nearest nodes only)->node assembly
  virtual void buildFaceElemToNodeGraph(const stk::mesh::PartVector&);// elem:face->node assembly
  virtual void buildNonConformalNodeGraph(const stk::mesh::PartVector&);// nonConformal->elem_node assembly
  virtual void buildOversetNodeGraph(const stk::mesh::PartVector&);// overset->elem_node assembly
  virtual void finalizeLinearSystem();

  /** Tag rows that must be handled as a Dirichlet BC node
   *
   *  @param[in] partVec List of parts that contain the Dirichlet nodes
   */
  virtual void buildDirichletNodeGraph(const stk::mesh::PartVector&);

  /** Tag rows that must be handled as a Dirichlet  node
   *
   *  @param[in] entities List of nodes where Dirichlet conditions are applied
   *
   *  \sa sierra::nalu::FixPressureAtNodeAlgorithm
   */
  virtual void buildDirichletNodeGraph(const std::vector<stk::mesh::Entity>&);

  // Matrix Assembly
  virtual void zeroSystem();

  virtual void sumInto(
      unsigned numEntities,
      const stk::mesh::Entity* entities,
      const SharedMemView<const double*> & rhs,
      const SharedMemView<const double**> & lhs,
      const SharedMemView<int*> & localIds,
      const SharedMemView<int*> & sortPermutation,
      const char * trace_tag
  );

  virtual void sumInto(
    const std::vector<stk::mesh::Entity> & sym_meshobj,
    std::vector<int> &scratchIds,
    std::vector<double> &scratchVals,
    const std::vector<double> & rhs,
    const std::vector<double> & lhs,
    const char *trace_tag
  );

  virtual void applyDirichletBCs(
    stk::mesh::FieldBase * solutionField,
    stk::mesh::FieldBase * bcValuesField,
    const stk::mesh::PartVector & parts,
    const unsigned beginPos,
    const unsigned endPos);

  /** Prepare assembly for overset fringe nodes
   *
   *  The overset fringe nodes are skipped over by the sumInto method during
   *  normal assembly process. This method toggles the flag to instruct sumInto
   *  that the constraint rows are being filled at this stage.
   */
  virtual void prepareConstraints(
    const unsigned,
    const unsigned)
  {
    checkSkippedRows_ = false;
  }

  /** Prepare assembly for Dirichlet-type rows
   *
   *  Dirichlet rows are skipped over by the sumInto method when the interior
   *  parts are processed. This method toggles the flag alerting the sumInto
   *  method that the Dirichlet rows will be processed next and sumInto can
   *  proceed.
   */
  virtual void resetRows(
    std::vector<stk::mesh::Entity>,
    const unsigned,
    const unsigned)
  {
    checkSkippedRows_ = false;
  }

  // Solve
  virtual int solve(stk::mesh::FieldBase * linearSolutionField);
  virtual void loadComplete();

  virtual void writeToFile(const char * filename, bool useOwned=true) {}
  virtual void writeSolutionToFile(const char * filename, bool useOwned=true) {}

protected:
  virtual void beginLinearSystemConstruction();

  int get_entity_hypre_id(const stk::mesh::Entity&);

  double copy_hypre_to_stk(stk::mesh::FieldBase*);

private:
  /** Flags indicating whether a particular row in the HYPRE matrix has been
   * filled or not.
   */
  enum RowFillStatus
  {
    RS_UNFILLED = 0, //!< Default status
    RS_FILLED        //!< sumInto filps to filled status once a row has been acted on
  };

  /** Flag indicating the type of row.
   *
   *  This flag is used to determine if the normal sumInto approach is used to
   *  populate the row, or a special method is used to handle that row. sumInto
   *  method will skip over the rows not marked RT_NORMAL and must be dealt with
   *  separately by other algorithms.
   */
  enum RowStatus
  {
    RT_NORMAL = 0, //!< A normal row that is summed into using sumInto
    RT_DIRICHLET,  //!< Rows with Dirichlet BC; no off-diagonal entries
    RT_OVERSET     //!< Overset fringe points; interpolation weights from other mesh
  };

  /** Dummy method to satisfy inheritance
   */
  void checkError(
    const int,
    const char*) {}

  //! The HYPRE matrix data structure
  mutable HYPRE_IJMatrix mat_;

  //! HYPRE right hand side data structure
  mutable HYPRE_IJVector rhs_;

  //! HYPRE solution vector
  mutable HYPRE_IJVector sln_;

  //! Track rows that have been updated during the assembly process
  std::vector<RowFillStatus> rowFilled_;

  //! Track the status of rows
  std::vector<RowStatus> rowStatus_;

  //! Track which rows are skipped
  std::unordered_set<int> skippedRows_;

  //! The lowest row owned by this MPI rank
  int iLower_;
  //! The highest row owned by this MPI rank
  int iUpper_;
  //! The lowest column owned by this MPI rank; currently jLower_ == iLower_
  int jLower_;
  //! The highest column owned by this MPI rank; currently jUpper_ == iUpper_
  int jUpper_;
  //! Total number of rows owned by this particular MPI rank
  int numRows_;
  //! Maximum Row ID in the Hypre linear system
  int maxRowID_;

  //! Flag indicating whether the linear system has been initialized
  bool systemInitialized_{false};

  //! Flag indicating whether the matrix has been assembled
  bool matrixFilled_{false};

  bool checkSkippedRows_{true};
};

}  // nalu
}  // sierra


#endif /* HYPRELINEARSYSTEM_H */
