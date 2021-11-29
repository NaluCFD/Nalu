/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNodalGradPAWElemAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include "SolutionOptions.h"
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleNodalGradPAWElemAlgorithm - balanced force Gjph interior
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradPAWElemAlgorithm::AssembleNodalGradPAWElemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *pressure,
  VectorFieldType *dpdx,
  const double buoyancyWeight)
  : Algorithm(realm, part),
    pressure_(pressure),
    density_(nullptr),
    interfaceCurvature_(NULL),
    surfaceTension_(NULL),
    vof_(NULL),
    dpdx_(dpdx),
    areaWeight_(nullptr),
    useShifted_(realm_.get_shifted_grad_op("pressure")),
    buoyancyWeight_(buoyancyWeight)
{
  // extract fields; nodal
  stk::mesh::MetaData & metaData = realm_.meta_data();
  density_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  interfaceCurvature_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "interface_curvature");
  surfaceTension_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "surface_tension");
  vof_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_of_fluid");
  areaWeight_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "png_area_weight");
  gravity_ = realm_.solutionOptions_->gravity_;
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradPAWElemAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // nodal fields to gather; gather everything other than what we are assembling
  std::vector<double> ws_pressure;
  std::vector<double> ws_density;  
  std::vector<double> ws_kappa;
  std::vector<double> ws_sigma;
  std::vector<double> ws_vof;

  std::vector<double> ws_coordinates;

  // geometry related to populate
  std::vector<double> ws_scs_areav;
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_det_j;
  std::vector<double> ws_shape_function;

  // integration point data that depends on size
  std::vector<double> dpdxIp(nDim);
  std::vector<double> dvofdxIp(nDim);
  double *p_dpdxIp = &dpdxIp[0];
  double *p_dvofdxIp = &dvofdxIp[0];

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_)  
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;
    const int *lrscv = meSCS->adjacentNodes();

    // algorithm related
    ws_pressure.resize(nodesPerElement);
    ws_density.resize(nodesPerElement);
    ws_kappa.resize(nodesPerElement);
    ws_sigma.resize(nodesPerElement);
    ws_vof.resize(nodesPerElement);
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_scs_areav.resize(numScsIp*nDim);
    ws_dndx.resize(nDim*numScsIp*nodesPerElement);
    ws_deriv.resize(nDim*numScsIp*nodesPerElement);
    ws_det_j.resize(numScsIp);
    ws_shape_function.resize(numScsIp*nodesPerElement);

    // pointers.
    double *p_pressure = &ws_pressure[0];
    double *p_density = &ws_density[0];
    double *p_kappa = &ws_kappa[0];
    double *p_sigma = &ws_sigma[0];
    double *p_vof = &ws_vof[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_scs_areav = &ws_scs_areav[0];
    double *p_dndx = &ws_dndx[0];
    double *p_shape_function = &ws_shape_function[0];

    meSCS->shape_fcn(&p_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        double * coords = stk::mesh::field_data(*coordinates, node);

        // gather scalars
        p_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
        p_density[ni] = *stk::mesh::field_data(*density_, node);
        p_kappa[ni]  = *stk::mesh::field_data(*interfaceCurvature_, node);
        p_sigma[ni]  = *stk::mesh::field_data(*surfaceTension_, node);
        p_vof[ni]  = *stk::mesh::field_data(*vof_, node);

        // gather vectors
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
        }
      }

      // compute geometry
      double scs_error = 0.0;
      meSCS->determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

      // compute dndx
      if (useShifted_)
        meSCS->shifted_grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);
      else
        meSCS->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);

      // start assembly
      for ( int ip = 0; ip < numScsIp; ++ip ) {

        // left and right nodes for this ip
        const int il = lrscv[2*ip];
        const int ir = lrscv[2*ip+1];

        stk::mesh::Entity nodeL = node_rels[il];
        stk::mesh::Entity nodeR = node_rels[ir];

        // pointer to fields to assemble
        double *gradPL = stk::mesh::field_data(*dpdx_, nodeL);
        double *gradPR = stk::mesh::field_data(*dpdx_, nodeR);
        double *areaWeightL = stk::mesh::field_data(*areaWeight_, nodeL);
        double *areaWeightR = stk::mesh::field_data(*areaWeight_, nodeR);

        // zero ip values
        double rhoIp = 0.0;
        double sigmaKappaIp = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          p_dpdxIp[j] = 0.0;
          p_dvofdxIp[j] = 0.0;
        }

        const int ipNpe = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function[ipNpe+ic];
          rhoIp += r*p_density[ic];
          sigmaKappaIp += r*p_sigma[ic]*p_kappa[ic];
          const double pressureIc = p_pressure[ic];
          const double vofIc = p_vof[ic];

          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            const double dxj = p_dndx[offSetDnDx+j];
            p_dpdxIp[j] += dxj*pressureIc;
            p_dvofdxIp[j] += dxj*vofIc;
          }
        }

        // assemble to il/ir - all addition here..
        const int ipNdim = ip*nDim;
        for ( int j = 0; j < nDim; ++j ) {
          const double absArea = std::abs(p_scs_areav[ipNdim+j]);
          double fac = absArea/rhoIp;
          gradPL[j] += fac*(p_dpdxIp[j] - rhoIp*gravity_[j]*buoyancyWeight_ - sigmaKappaIp*dvofdxIp[j]);
          gradPR[j] += fac*(p_dpdxIp[j] - rhoIp*gravity_[j]*buoyancyWeight_ - sigmaKappaIp*dvofdxIp[j]);
          areaWeightL[j] += absArea;
          areaWeightR[j] += absArea;
        }
      }
    }
  }
}
  
} // namespace nalu
} // namespace Sierra
