/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleScalarEdgeDiffContactSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>

#include <ContactInfo.h>
#include <ContactManager.h>
#include <HermitePolynomialInterpolation.h>
#include <FieldTypeDef.h>
#include <HaloInfo.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <TimeIntegrator.h>

#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleScalarEdgeDiffContactSolverAlgorithm - add LHS/RHS for q
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEdgeDiffContactSolverAlgorithm::AssembleScalarEdgeDiffContactSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  VectorFieldType *dqdx,
  ScalarFieldType *diffFluxCoeff,
  bool useHermiteInterpolation)
  : SolverAlgorithm(realm, part, eqSystem),
    hermite_(NULL),
    scalarQ_(scalarQ),
    dqdx_(dqdx),
    diffFluxCoeff_(diffFluxCoeff)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  haloMdot_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_mdot");

  // populate fieldVec
  ghostFieldVec_.push_back(scalarQ_);
  ghostFieldVec_.push_back(dqdx_);
  ghostFieldVec_.push_back(diffFluxCoeff_);
  ghostFieldVec_.push_back(coordinates_);

  if ( useHermiteInterpolation ) {
    const int nDim = meta_data.spatial_dimension();
    if ( nDim == 2 ) {
      hermite_ = new HermitePolynomialInterpolationFourPoint();
    }
    else {
      hermite_ = new HermitePolynomialInterpolationEightPoint();
    }
  }
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEdgeDiffContactSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildEdgeHaloNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEdgeDiffContactSolverAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // space for LHS/RHS (nodesPerElem+1)*(nodesPerElem+1); nodesPerElem+1
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // deal with state
  ScalarFieldType &scalarQNp1 = scalarQ_->field_of_state(stk::mesh::StateNP1);

  // space for interpolated right state (halo)
  double qNp1R;
  double diffFluxCoeffR;
  std::vector<double> dqdxR(nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
  const int sizeOfVectorField = nDim;
  
  // parallel communicate ghosted entities
  if ( NULL != realm_.contactManager_->contactGhosting_ )
    stk::mesh::communicate_field_data(*(realm_.contactManager_->contactGhosting_), ghostFieldVec_);
  
  // iterate contactInfoVec_
  std::vector<ContactInfo *>::iterator ii;
  for( ii=realm_.contactManager_->contactInfoVec_.begin();
       ii!=realm_.contactManager_->contactInfoVec_.end(); ++ii ) {
    
    // get master element type for this contactInfo
    MasterElement *meSCS  = (*ii)->meSCS_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    std::vector <double > elemNodalQ(nodesPerElement);
    std::vector <double > elemNodalDiffFluxCoeff(nodesPerElement);
    std::vector <double > elemNodalCoords(nDim*nodesPerElement);
    std::vector <double > elemNodalDqdx(nDim*nodesPerElement);
    std::vector <double > shpfc(nodesPerElement);

    // resize some things; matrix related
    const int npePlusOne = nodesPerElement+1;
    const int lhsSize = npePlusOne*npePlusOne;
    const int rhsSize = npePlusOne;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(npePlusOne);

    // pointer to lhs/rhs
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];

    // iterate halo face nodes
    std::map<uint64_t, HaloInfo *>::iterator iterHalo;
    for (iterHalo  = (*ii)->haloInfoMap_.begin();
         iterHalo != (*ii)->haloInfoMap_.end();
         ++iterHalo) {

      // halo info object of interest
      HaloInfo * infoObject = (*iterHalo).second;

      // zeroing of lhs/rhs
      for ( int k = 0; k < lhsSize; ++k ) {
        p_lhs[k] = 0.0;
      }
      for ( int k = 0; k < rhsSize; ++k ) {
        p_rhs[k] = 0.0;
      }

      // pointer to edge area vector
      const double *p_areaVec = &infoObject->haloEdgeAreaVec_[0];

      // extract element mesh object and global id for face node
      stk::mesh::Entity elem  = infoObject->owningElement_;

      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(elem);
      const int num_nodes = bulk_data.num_nodes(elem);

      // now load the elemental values for future interpolation; fill in connected nodes
      connected_nodes[0] = infoObject->faceNode_;
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        connected_nodes[ni+1] = node;

        elemNodalQ[ni] = *stk::mesh::field_data(*scalarQ_, node);
        elemNodalDiffFluxCoeff[ni] = *stk::mesh::field_data(*diffFluxCoeff_, node);

        // load up vectors
        const double * Gjq = stk::mesh::field_data(*dqdx_, node );
        const double * coords = stk::mesh::field_data(*coordinates_, node );
        for ( int j = 0; j < nDim; ++j ) {
          const int offSet = j*nodesPerElement +ni;
          elemNodalDqdx[offSet] = Gjq[j];
          elemNodalCoords[offSet] = coords[j];
        }
      }

      // extract nodal fields; right state is Halo and requires inperpolation
      const double *coordL = stk::mesh::field_data(*coordinates_, infoObject->faceNode_);
      const double *coordR = &infoObject->haloNodalCoords_[0];

      const double qNp1L = *stk::mesh::field_data(scalarQNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalQ[0],
        &qNp1R);

      // possible Hermite polynomial interpolation
      if (NULL != hermite_) {
        double qNp1H = 0;
        hermite_->do_hermite(&elemNodalQ[0], &elemNodalCoords[0], &elemNodalDqdx[0], &coordR[0], qNp1H);
        qNp1R = qNp1H;
      }

      const double diffFluxCoeffL = *stk::mesh::field_data(*diffFluxCoeff_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalDiffFluxCoeff[0],
        &diffFluxCoeffR);

      // deal with nodal dqdx
      const double *dqdxL = stk::mesh::field_data(*dqdx_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfVectorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalDqdx[0],
        &(dqdxR[0]));

      // ip props
      const double viscIp = 0.5*(diffFluxCoeffL + diffFluxCoeffR);

      // compute geometry
      double axdx = 0.0;
      double asq = 0.0;
      for (int j = 0; j < nDim; ++j ) {
        const double dxj = coordR[j] - coordL[j];
        const double axj = p_areaVec[j];
        axdx += axj*dxj;
        asq += axj*axj;
      }

      const double inv_axdx = 1.0/axdx;

      // NOC
      double nonOrth = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        const double axj = p_areaVec[j];
        const double dxj = coordR[j] - coordL[j];
        // now non-orth (over-relaxed procedure of Jasek)
        const double kxj = axj - asq*inv_axdx*dxj;
        const double GjIp = 0.5*(dqdxL[j] + dqdxR[j]);
        nonOrth += -viscIp*kxj*GjIp;
      }

      // iterate owning element nodes and form pair
      meSCS->general_shape_fcn(1, &(infoObject->isoParCoords_[0]), &shpfc[0]);

      const double lhsFac = -viscIp*asq*inv_axdx;
      const double diffFlux = lhsFac*(qNp1R-qNp1L) + nonOrth;

      // setup for LHS; row left is easy - simply zero

      // left node (face)
      p_rhs[0] = -diffFlux;
      p_lhs[0] = -lhsFac;

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        p_lhs[ni+1] = lhsFac*shpfc[ni];
      }

      // apply to linear system
      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
    }
  }
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEdgeDiffContactSolverAlgorithm::~AssembleScalarEdgeDiffContactSolverAlgorithm()
{
  if ( NULL != hermite_ )
    delete hermite_;
}

} // namespace nalu
} // namespace Sierra
