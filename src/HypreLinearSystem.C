/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifdef NALU_USES_HYPRE

#include "HypreLinearSystem.h"
#include "HypreDirectSolver.h"
#include "Realm.h"
#include "EquationSystem.h"
#include "LinearSolver.h"
#include "PeriodicManager.h"
#include "NonConformalManager.h"
#include "overset/OversetManager.h"
#include "overset/OversetInfo.h"


#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"

#include <cmath>

namespace sierra {
namespace nalu {

HypreLinearSystem::HypreLinearSystem(
  Realm& realm,
  const unsigned numDof,
  EquationSystem* eqSys,
  LinearSolver* linearSolver)
  : LinearSystem(realm, numDof, eqSys, linearSolver),
    rowFilled_(0),
    rowStatus_(0)
{
}

HypreLinearSystem::~HypreLinearSystem()
{
  // if (linearSolver_ != nullptr)
  //   linearSolver_->destroyLinearSolver();

  if (systemInitialized_) {
    HYPRE_IJMatrixDestroy(mat_);
    HYPRE_IJVectorDestroy(rhs_);
    HYPRE_IJVectorDestroy(sln_);
    systemInitialized_ = false;
  }
}

void
HypreLinearSystem::beginLinearSystemConstruction()
{
  if (inConstruction_) return;
  inConstruction_ = true;

  const int rank = realm_.bulk_data().parallel_rank();

  if (rank == 0) {
    iLower_ = realm_.hypreILower_;
  } else {
    iLower_ = realm_.hypreILower_ * numDof_ ;
  }

  iUpper_ = realm_.hypreIUpper_  * numDof_ - 1;
  // For now set column indices the same as row indices
  jLower_ = iLower_;
  jUpper_ = iUpper_;

  // The total number of rows handled by this MPI rank for Hypre
  numRows_ = (iUpper_ - iLower_ + 1);
  // Total number of global rows in the system
  maxRowID_ = realm_.hypreNumNodes_ * numDof_ - 1;

#if 0
  if (numDof_ > 0)
    std::cerr << rank << "\t" << numDof_ << "\t"
              << realm_.hypreILower_ << "\t" << realm_.hypreIUpper_ << "\t"
                << iLower_ << "\t" << iUpper_ << "\t"
                << numRows_ << "\t" << maxRowID_ << std::endl;
#endif
  // Allocate memory for the arrays used to track row types and row filled status.
  rowFilled_.resize(numRows_);
  rowStatus_.resize(numRows_);
  skippedRows_.clear();
  // All nodes start out as NORMAL; "build*NodeGraph" methods might alter the
  // row status to modify behavior of sumInto method.
  for (int i=0; i < numRows_; i++)
    rowStatus_[i] = RT_NORMAL;

  auto& bulk = realm_.bulk_data();
  std::vector<const stk::mesh::FieldBase*> fVec{realm_.hypreGlobalId_};

  stk::mesh::copy_owned_to_shared(bulk, fVec);
  stk::mesh::communicate_field_data(bulk.aura_ghosting(), fVec);

  if (realm_.oversetManager_ != nullptr &&
      realm_.oversetManager_->oversetGhosting_ != nullptr)
    stk::mesh::communicate_field_data(
      *realm_.oversetManager_->oversetGhosting_, fVec);

  if (realm_.nonConformalManager_ != nullptr &&
      realm_.nonConformalManager_->nonConformalGhosting_ != nullptr)
    stk::mesh::communicate_field_data(
      *realm_.nonConformalManager_->nonConformalGhosting_, fVec);

  if (realm_.periodicManager_ != nullptr &&
      realm_.periodicManager_->periodicGhosting_ != nullptr) {
    realm_.periodicManager_->parallel_communicate_field(realm_.hypreGlobalId_);
    realm_.periodicManager_->periodic_parallel_communicate_field(
      realm_.hypreGlobalId_);
  }
}

void
HypreLinearSystem::buildNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();
}

void
HypreLinearSystem::buildFaceToNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();
}

void
HypreLinearSystem::buildEdgeToNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();
}

void
HypreLinearSystem::buildElemToNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();
}

void
HypreLinearSystem::buildFaceElemToNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();
}

void
HypreLinearSystem::buildReducedElemToNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();
}

void
HypreLinearSystem::buildNonConformalNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();
}

void
HypreLinearSystem::buildOversetNodeGraph(
  const stk::mesh::PartVector&)
{
  beginLinearSystemConstruction();

  // Mark all the fringe nodes as skipped so that sumInto doesn't add into these
  // rows during assembly process
  for(auto* oinfo: realm_.oversetManager_->oversetInfoVec_) {
    auto node = oinfo->orphanNode_;
    int hid = *stk::mesh::field_data(*realm_.hypreGlobalId_, node);
    for (size_t d=0; d < numDof_; d++) {
      int lid = hid * numDof_ + d;
      // rowStatus_[lid] = RT_OVERSET;
      skippedRows_.insert(lid);
    }
  }
}

void
HypreLinearSystem::buildDirichletNodeGraph(
  const stk::mesh::PartVector& parts)
{
  beginLinearSystemConstruction();

  // Grab nodes regardless of whether they are owned or shared
  const stk::mesh::Selector sel = stk::mesh::selectUnion(parts);
  const auto& bkts = realm_.get_buckets(
    stk::topology::NODE_RANK, sel);

  for (auto b: bkts) {
    for (size_t in=0; in < b->size(); in++) {
      auto node = (*b)[in];
      int hid = *stk::mesh::field_data(*realm_.hypreGlobalId_, node);
      for (size_t d=0; d < numDof_; d++) {
        int lid = hid * numDof_ + d;
        // rowStatus_[lid] = RT_DIRICHLET;
        skippedRows_.insert(lid);
      }
    }
  }
}

void
HypreLinearSystem::buildDirichletNodeGraph(
  const std::vector<stk::mesh::Entity>& nodeList)
{
  beginLinearSystemConstruction();

  for (const auto& node: nodeList) {
    int hid = get_entity_hypre_id(node);
    for (size_t d=0; d < numDof_; d++) {
      int lid = hid * numDof_ + d;
      // rowStatus_[lid] = RT_DIRICHLET;
      skippedRows_.insert(lid);
    }
  }
}

void
HypreLinearSystem::finalizeLinearSystem()
{
  ThrowRequire(inConstruction_);
  inConstruction_ = false;

  // Prepare for matrix assembly and set all entry flags to "unfilled"
  for (int i=0; i < numRows_; i++)
    rowFilled_[i] = RS_UNFILLED;

  MPI_Comm comm = realm_.bulk_data().parallel();
  // Now perform HYPRE assembly so that the data structures are ready to be used
  // by the solvers/preconditioners.
  HypreDirectSolver* solver = reinterpret_cast<HypreDirectSolver*>(linearSolver_);

  HYPRE_IJMatrixCreate(comm, iLower_, iUpper_, jLower_, jUpper_, &mat_);
  HYPRE_IJMatrixSetObjectType(mat_, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(mat_);
  HYPRE_IJMatrixGetObject(mat_, (void**)&(solver->parMat_));

  HYPRE_IJVectorCreate(comm, iLower_, iUpper_, &rhs_);
  HYPRE_IJVectorSetObjectType(rhs_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(rhs_);
  HYPRE_IJVectorGetObject(rhs_, (void**)&(solver->parRhs_));

  HYPRE_IJVectorCreate(comm, iLower_, iUpper_, &sln_);
  HYPRE_IJVectorSetObjectType(sln_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(sln_);
  HYPRE_IJVectorGetObject(sln_, (void**)&(solver->parSln_));

  // At this stage the LHS and RHS data structures are ready for
  // sumInto/assembly.
  systemInitialized_ = true;
}

void
HypreLinearSystem::loadComplete()
{
  // All algorithms have called sumInto and populated LHS/RHS. Now we are ready
  // to finalize the matrix at the HYPRE end. However, before we do that we need
  // to process unfilled rows and process them appropriately. Any row acted on
  // by sumInto method will have toggled the rowFilled_ array to RS_FILLED
  // status. Before finalizing assembly, we process rows that still have an
  // RS_UNFILLED status and set their diagonal entries to 1.0 (dummy row)
  //
  // TODO: Alternate design to eliminate dummy rows. This will require
  // load-balancing on HYPRE end.

  int hnrows = 1;
  int hncols = 1;
  double getval;
  double setval = 1.0;
  for (int i=0; i < numRows_; i++) {
    if (rowFilled_[i] == RS_FILLED) continue;
    int lid = iLower_ + i;
    HYPRE_IJMatrixGetValues(mat_, hnrows, &hncols, &lid, &lid, &getval);
    if (std::fabs(getval) < 1.0e-12)
      HYPRE_IJMatrixSetValues(mat_, hnrows, &hncols, &lid, &lid, &setval);
  }

  // Now perform HYPRE assembly so that the data structures are ready to be used
  // by the solvers/preconditioners.
  HypreDirectSolver* solver = reinterpret_cast<HypreDirectSolver*>(linearSolver_);

  HYPRE_IJMatrixAssemble(mat_);
  HYPRE_IJMatrixGetObject(mat_, (void**)&(solver->parMat_));

  HYPRE_IJVectorAssemble(rhs_);
  HYPRE_IJVectorGetObject(rhs_, (void**)&(solver->parRhs_));

  HYPRE_IJVectorAssemble(sln_);
  HYPRE_IJVectorGetObject(sln_, (void**)&(solver->parSln_));

  solver->comm_ = realm_.bulk_data().parallel();
}

void
HypreLinearSystem::zeroSystem()
{
  MPI_Comm comm = realm_.bulk_data().parallel();
  HypreDirectSolver* solver = reinterpret_cast<HypreDirectSolver*>(linearSolver_);

  if (systemInitialized_)
    HYPRE_IJMatrixDestroy(mat_);
  HYPRE_IJMatrixCreate(comm, iLower_, iUpper_, jLower_, jUpper_, &mat_);
  HYPRE_IJMatrixSetObjectType(mat_, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(mat_);
  HYPRE_IJMatrixGetObject(mat_, (void**)&(solver->parMat_));

  HYPRE_IJVectorInitialize(rhs_);
  HYPRE_IJVectorInitialize(sln_);

  HYPRE_ParVectorSetConstantValues(solver->parRhs_, 0.0);
  HYPRE_ParVectorSetConstantValues(solver->parSln_, 0.0);

  // Prepare for matrix assembly and set all entry flags to "unfilled"
  for (int i=0; i < numRows_; i++)
    rowFilled_[i] = RS_UNFILLED;

  // Reset overset flag so that sumInto only processes non-overset fringe rows
  // until we are ready to process overset constraint rows
  checkSkippedRows_ = true;
}

void
HypreLinearSystem::sumInto(
  unsigned numEntities,
  const stk::mesh::Entity* entities,
  const SharedMemView<const double*>& rhs,
  const SharedMemView<const double**>& lhs,
  const SharedMemView<int*>& localIds,
  const SharedMemView<int*>& sortPermutations,
  const char* trace_tag)
{
  const size_t n_obj = numEntities;
  int numRows = n_obj * numDof_;

  ThrowAssertMsg(lhs.is_contiguous(), "LHS assumed contiguous");
  ThrowAssertMsg(rhs.is_contiguous(), "RHS assumed contiguous");
  ThrowAssertMsg(localIds.is_contiguous(), "localIds assumed contiguous");
  ThrowAssertMsg(sortPermutation.is_contiguous(), "sortPermutation assumed contiguous");

  for (size_t in=0; in < n_obj; in++) {
    int hid = get_entity_hypre_id(entities[in]);
    int localOffset = hid * numDof_;
    for (size_t d=0; d < numDof_; d++) {
      size_t lid = in * numDof_ + d;
      localIds[lid] = localOffset + d;
    }
  }

  for (int ir=0; ir < numRows; ir++) {
      int lid = localIds[ir];
      auto it = skippedRows_.find(lid);
      if (checkSkippedRows_ && (it != skippedRows_.end())) continue;

      const double* cur_lhs = &lhs(ir, 0);
      HYPRE_IJMatrixAddToValues(mat_, 1, &numRows, &lid,
                                &localIds[0], cur_lhs);
      HYPRE_IJVectorAddToValues(rhs_, 1, &lid, &rhs[ir]);

      if ((lid >= iLower_) && (lid <= iUpper_))
        rowFilled_[lid - iLower_] = RS_FILLED;
  }
}

void
HypreLinearSystem::sumInto(
  const std::vector<stk::mesh::Entity>& entities,
  std::vector<int>& scratchIds,
  std::vector<double>& scratchVals,
  const std::vector<double>& rhs,
  const std::vector<double>& lhs,
  const char* trace_tag)
{
  const size_t n_obj = entities.size();
  int numRows = n_obj * numDof_;

  ThrowAssert(numRows == rhs.size());
  ThrowAssert(numRows*numRows == lhs.size());

  for (size_t in=0; in < n_obj; in++) {
    int hid = get_entity_hypre_id(entities[in]);
    int localOffset = hid * numDof_;
    for (size_t d=0; d < numDof_; d++) {
      size_t lid = in * numDof_ + d;
      scratchIds[lid] = localOffset + d;
    }
  }

  for (int ir=0; ir < numRows; ir++) {
    int lid = scratchIds[ir];
    auto it = skippedRows_.find(lid);
    if (checkSkippedRows_ && (it != skippedRows_.end())) continue;

    for (int c=0; c < numRows; c++)
      scratchVals[c] = lhs[ir * numRows + c];

    HYPRE_IJMatrixAddToValues(mat_, 1, &numRows, &lid,
                              &scratchIds[0], &scratchVals[0]);
    HYPRE_IJVectorAddToValues(rhs_, 1, &lid, &rhs[ir]);
    if ((lid >= iLower_) && (lid <= iUpper_))
        rowFilled_[lid - iLower_] = RS_FILLED;
  }
}

void
HypreLinearSystem::applyDirichletBCs(
  stk::mesh::FieldBase* solutionField,
  stk::mesh::FieldBase* bcValuesField,
  const stk::mesh::PartVector& parts,
  const unsigned,
  const unsigned)
{
  auto& meta = realm_.meta_data();

  const stk::mesh::Selector sel = (
    meta.locally_owned_part() &
    stk::mesh::selectUnion(parts) &
    stk::mesh::selectField(*solutionField) &
    !(realm_.get_inactive_selector()));

  const auto& bkts = realm_.get_buckets(
    stk::topology::NODE_RANK, sel);

  int ncols = 1;
  double diag_value = 1.0;
  for (auto b: bkts) {
    const double* solution = (double*)stk::mesh::field_data(
      *solutionField, *b);
    const double* bcValues = (double*)stk::mesh::field_data(
      *bcValuesField, *b);

    for (size_t in=0; in < b->size(); in++) {
      auto node = (*b)[in];
      int hid = *stk::mesh::field_data(*realm_.hypreGlobalId_, node);

      for (size_t d=0; d<numDof_; d++) {
        int lid = hid * numDof_ + d;
        double bcval = bcValues[in*numDof_ + d] - solution[in*numDof_ + d];

        HYPRE_IJMatrixSetValues(mat_, 1, &ncols, &lid, &lid, &diag_value);
        HYPRE_IJVectorSetValues(rhs_, 1, &lid, &bcval);
        rowFilled_[lid - iLower_] = RS_FILLED;
      }
    }
  }
}

int
HypreLinearSystem::get_entity_hypre_id(const stk::mesh::Entity& node)
{
  auto& bulk = realm_.bulk_data();
  const auto naluId = *stk::mesh::field_data(*realm_.naluGlobalId_, node);
  const auto mnode = bulk.get_entity(stk::topology::NODE_RANK, naluId);
  if (!bulk.is_valid(node))
    throw std::runtime_error("BAD STK NODE");
  int hid = *stk::mesh::field_data(*realm_.hypreGlobalId_, mnode);

  if ((hid < 0) || (hid > maxRowID_)) {
    std::cerr << bulk.parallel_rank() << "\t"
              << hid << "\t" << iLower_ << "\t" << iUpper_ << std::endl;
    throw std::runtime_error("BAD STK to hypre conversion");
  }

  return hid;
}

int
HypreLinearSystem::solve(stk::mesh::FieldBase* linearSolutionField)
{
  HypreDirectSolver* solver = reinterpret_cast<HypreDirectSolver*>(
    linearSolver_);

  if (solver->getConfig()->getWriteMatrixFiles()) {
    const std::string matFile = eqSysName_ + ".IJM.mat.";
    const std::string rhsFile = eqSysName_ + ".IJV.rhs.";
    HYPRE_IJMatrixPrint(mat_, matFile.c_str());
    HYPRE_IJVectorPrint(rhs_, rhsFile.c_str());
  }

  int iters = 0;
  double finalResidNorm = 0.0;

  // Call solve
  int status = 0;

  status = solver->solve(iters, finalResidNorm);

  if (solver->getConfig()->getWriteMatrixFiles()) {
    const std::string slnFile = eqSysName_ + ".IJV.sln.";
    HYPRE_IJVectorPrint(sln_, slnFile.c_str());
  }

  double norm2 = copy_hypre_to_stk(linearSolutionField);
  sync_field(linearSolutionField);

  linearSolveIterations_ = iters;
  linearResidual_ = finalResidNorm;
  nonLinearResidual_ = realm_.l2Scaling_ * norm2;

  if (eqSys_->firstTimeStepSolve_)
    firstNonLinearResidual_ = nonLinearResidual_;

  scaledNonLinearResidual_ =
    nonLinearResidual_ /
    std::max(std::numeric_limits<double>::epsilon(), firstNonLinearResidual_);

  if (provideOutput_) {
    const int nameOffset = eqSysName_.length() + 8;
    NaluEnv::self().naluOutputP0()
      << std::setw(nameOffset) << std::right << eqSysName_
      << std::setw(32 - nameOffset) << std::right << iters << std::setw(18)
      << std::right << finalResidNorm << std::setw(15) << std::right
      << nonLinearResidual_ << std::setw(14) << std::right
      << scaledNonLinearResidual_ << std::endl;
  }

  eqSys_->firstTimeStepSolve_ = false;

  return status;
}

double
HypreLinearSystem::copy_hypre_to_stk(
  stk::mesh::FieldBase* stkField)
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  const auto sel = stk::mesh::selectField(*stkField)
    & meta.locally_owned_part()
    & !(stk::mesh::selectUnion(realm_.get_slave_part_vector()))
    & !(realm_.get_inactive_selector());

  const auto& bkts = bulk.get_buckets(
    stk::topology::NODE_RANK, sel);

  double lclnorm2 = 0.0;
  for (auto b: bkts) {
    double* field = (double*) stk::mesh::field_data(*stkField, *b);
    for (size_t in=0; in < b->size(); in++) {
      auto node = (*b)[in];
      int hid = get_entity_hypre_id(node);

      for (size_t d=0; d < numDof_; d++) {
        int lid = hid * numDof_ + d;
        int sid = in * numDof_ + d;
        HYPRE_IJVectorGetValues(sln_, 1, &lid, &field[sid]);

        lclnorm2 += field[sid] * field[sid];
      }
    }
  }

  double gblnorm2 = 0.0;
  stk::all_reduce_sum(bulk.parallel(), &lclnorm2, &gblnorm2, 1);

  return std::sqrt(gblnorm2);
}

}  // nalu
}  // sierra

#endif
