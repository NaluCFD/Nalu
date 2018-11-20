/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include "ComputeDynamicPressureAlgorithm.h"
#include "Algorithm.h"

#include "FieldTypeDef.h"
#include "Realm.h"
#include "SolutionOptions.h"
#include "master_element/MasterElement.h"

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ComputeDynamicPressureAlgorithm - mdot continuity open bc
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeDynamicPressureAlgorithm::ComputeDynamicPressureAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  const bool useShifted)
  : Algorithm(realm, part),
    useShifted_(useShifted),
    density_(NULL),
    openMassFlowRate_(NULL),
    exposedAreaVec_(NULL),
    dynamicPressure_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  openMassFlowRate_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "open_mass_flow_rate");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  dynamicPressure_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "dynamic_pressure");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeDynamicPressureAlgorithm::~ComputeDynamicPressureAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeDynamicPressureAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  
  // nodal fields to gather
  std::vector<double> ws_density;

  // master element
  std::vector<double> ws_face_shape_function;

  // deal with state
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = b.topology().num_nodes();
    const int numScsBip = meFC->numIntPoints_;

    // algorithm related; element
    ws_density.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_density = &ws_density[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape functions; boundary
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < b.size(); ++k ) {
      
      // get face
      stk::mesh::Entity face = b[k];

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      int num_face_nodes = bulk_data.num_nodes(face);
      // sanity check on num nodes
      ThrowAssert( num_face_nodes == nodesPerFace );
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];

        // gather scalars
        p_density[ni]    = *stk::mesh::field_data(densityNp1, node);
      }

      // pointer to face data; will scatter dynamicP
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      const double * mdot = stk::mesh::field_data(*openMassFlowRate_, face);
      double * dynamicP = stk::mesh::field_data(*dynamicPressure_, face);

      // loop over boundary ips
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int ipNdim = ip*nDim;

        // compute area magnitude
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[ipNdim+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // interpolate density to bip
        double rhoBip = 0.0;
        const int ipNpf = ip*nodesPerFace;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[ipNpf+ic];
          rhoBip += r*p_density[ic];
        }
        
        // compute dynamic pressure
        const double uBip = mdot[ip]/(rhoBip*aMag);
        dynamicP[ip] = uBip < 0.0 ? 0.5*rhoBip*uBip*uBip : 0.0;
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
