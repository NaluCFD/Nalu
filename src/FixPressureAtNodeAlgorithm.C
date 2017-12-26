/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "FixPressureAtNodeAlgorithm.h"
#include "FixPressureAtNodeInfo.h"
#include "Realm.h"
#include "EquationSystem.h"
#include "TpetraLinearSystem.h"
#include "SolutionOptions.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <limits>

namespace sierra {
namespace nalu {

FixPressureAtNodeAlgorithm::FixPressureAtNodeAlgorithm(
  Realm& realm,
  stk::mesh::Part* part,
  EquationSystem* eqSystem)
  : SolverAlgorithm(realm, part, eqSystem),
    info_(*(realm.solutionOptions_->fixPressureInfo_)),
    meshMotion_(realm.does_mesh_move())
{
  auto& meta = realm_.meta_data();

  coordinates_ = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, realm_.get_coordinates_name());
  pressure_ = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "pressure");
}


FixPressureAtNodeAlgorithm::~FixPressureAtNodeAlgorithm()
{}

void
FixPressureAtNodeAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildDirichletNodeGraph(refNodeList_);
}

void
FixPressureAtNodeAlgorithm::execute()
{
  if (doInit_)
    initialize();

  // Reset LHS and RHS for this matrix
  eqSystem_->linsys_->resetRows(refNodeList_, 0, 1);

  int numNodes = refNodeList_.size();
  ThrowAssertMsg(numNodes <= 1,
                 "Invalid number of nodes encountered in FixPressureAtNodeAlgorithm");

  // Fix the pressure for this node only if this is proc is owner
  if (numNodes > 0 && fixPressureNode_) {
    std::vector<double> lhs(1);
    std::vector<double> rhs(1);
    std::vector<int> scratchIds(1);
    std::vector<double> scratchVals(1);
    stk::mesh::Entity node = refNodeList_[0];
    const double pressureN = *stk::mesh::field_data(*pressure_, node);

    lhs[0] = 1.0; // Set diagonal entry to 1.0
    rhs[0] = info_.refPressure_ - pressureN;

    apply_coeff(refNodeList_, scratchIds, scratchVals, rhs, lhs, __FILE__);
  }
}

void
FixPressureAtNodeAlgorithm::initialize()
{
  // Reset data structures
  refNodeList_.clear();

  if (info_.lookupType_ == FixPressureAtNodeInfo::SPATIAL_LOCATION) {
    // Determine the nearest node where pressure is referenced
    auto nodeID = determine_nearest_node();
    process_pressure_fix_node(nodeID);
  }
  else if (info_.lookupType_ == FixPressureAtNodeInfo::STK_NODE_ID) {
    process_pressure_fix_node(info_.stkNodeId_);
  }

  // Flip init flag
  doInit_ = false;
}

stk::mesh::EntityId
FixPressureAtNodeAlgorithm::determine_nearest_node()
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  int nDim = meta.spatial_dimension();
  auto& refLoc = info_.location_;

  // Get the target search part vector
  auto& partNames = info_.searchParts_;
  auto nParts = partNames.size();
  stk::mesh::PartVector parts(nParts);
  for (size_t i=0; i<nParts; i++) {
    stk::mesh::Part* part = meta.get_part(partNames[i]);
    if ( nullptr != part)
      parts[i] = part;
    else
      throw std::runtime_error("FixPressureAtNodeAlgorithm: Target search part is null " +
                               partNames[i]);
  }

  // Determine the nearest node in this processor
  stk::mesh::Entity nearestNode;
  stk::mesh::Selector sel = meta.locally_owned_part() &
    stk::mesh::selectUnion(parts);
  auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, sel);

  double distSqr = std::numeric_limits<double>::max();
  for (auto b: buckets) {
    auto length = b->size();

    for (size_t i=0; i < length; i++) {
      auto node = (*b)[i];
      const double* coords = stk::mesh::field_data(*coordinates_, node);

      double dist = 0.0;
      for (int j=0; j<nDim; j++) {
        double xdiff = (refLoc[j]-coords[j]);
        dist += xdiff * xdiff;
      }
      if (dist < distSqr) {
        distSqr = dist;
        nearestNode = node;
      }
    }
  }

  // Determine the global minimum
  std::vector<double> minDistList(bulk.parallel_size());
  MPI_Allgather(&distSqr, 1, MPI_DOUBLE, minDistList.data(), 1, MPI_DOUBLE,
                bulk.parallel());
  int minDistProc = -1;
  double minDist = std::numeric_limits<double>::max();
  for (int i=0; i < bulk.parallel_size(); i++) {
    if (minDistList[i] < minDist) {
      minDist = minDistList[i];
      minDistProc = i;
    }
  }

  // Communicate the nearest node ID to all processors.
  stk::mesh::EntityId nodeID = 0;
  stk::mesh::EntityId g_nodeID;
  if (minDistProc == bulk.parallel_rank())
    nodeID = bulk.identifier(nearestNode);
  stk::all_reduce_max(bulk.parallel(), &nodeID, &g_nodeID, 1);

  return g_nodeID;
}

void
FixPressureAtNodeAlgorithm::process_pressure_fix_node(
  const stk::mesh::EntityId nodeID)
{
  auto& bulk = realm_.bulk_data();

  // Store the target node on the owning processor as well as the shared processors.
  stk::mesh::Entity targetNode = bulk.get_entity(stk::topology::NODE_RANK, nodeID);
  if (bulk.is_valid(targetNode) &&
      (bulk.bucket(targetNode).owned() ||
       bulk.bucket(targetNode).shared())) {
    refNodeList_.push_back(targetNode);

    // Only apply pressure correction on the owning processor
    fixPressureNode_ = bulk.bucket(targetNode).owned();

#if 0
    if (bulk.parallel_owner_rank(targetNode) == bulk.parallel_rank()) {
      std::cerr
        << "FixPressureAtNodeAlgorithm: Node ID = "
        << nodeID << "; Proc ID = " << bulk.parallel_rank() << std::endl;
    }
#endif
  } else {
    fixPressureNode_ = false;
  }
}


}  // nalu
}  // sierra
