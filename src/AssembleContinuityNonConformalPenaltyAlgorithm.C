/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleContinuityNonConformalPenaltyAlgorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <cmath>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleContinuityNonConformalPenaltyAlgorithm - nodal lambda, flux
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleContinuityNonConformalPenaltyAlgorithm::AssembleContinuityNonConformalPenaltyAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *ncNormalFlux,
  ScalarFieldType *ncPenalty,
  ScalarFieldType *ncArea,
  const bool useShifted)
  : Algorithm(realm, part),
    ncNormalFlux_(ncNormalFlux),
    ncPenalty_(ncPenalty),
    ncArea_(ncArea),
    useShifted_(useShifted),
    meshMotion_(realm_.does_mesh_move()),
    velocityRTM_(NULL),
    density_(NULL),
    coordinates_(NULL),
    exposedAreaVec_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
 if ( meshMotion_ )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleContinuityNonConformalPenaltyAlgorithm::execute()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;

  // deal with interpolation procedure
  const double interpTogether = realm_.get_mdot_interp();
  const double om_interpTogether = 1.0-interpTogether;

  // fixed size
  std::vector<double> uBip(nDim);
  std::vector<double> rho_uBip(nDim);

  // pointers to fixed values
  double *p_uBip = &uBip[0];
  double *p_rho_uBip = &rho_uBip[0];
 
  // nodal fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_density;
  std::vector<double> ws_vrtm;

  // master element
  std::vector<double> ws_face_shape_function;

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // extract connected element topology
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    ThrowAssert ( parentTopo.size() == 1 );
    stk::topology theElemTopo = parentTopo[0];
    MasterElement *meSCS = realm_.get_surface_master_element(theElemTopo);
    const int nodesPerElement = meSCS->nodesPerElement_;

    // face master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;

    // size some things that are useful
    int num_face_nodes = b.topology().num_nodes();
    std::vector<int> face_node_ordinals(num_face_nodes);

    // algorithm related; element nodes
    ws_coordinates.resize(nodesPerElement*nDim);
   
    // face nodes
    ws_density.resize(nodesPerFace);
    ws_vrtm.resize(nodesPerElement*nDim);
    ws_face_shape_function.resize(nodesPerFace*nodesPerFace);
  
    // pointers
    double *p_coordinates = &ws_coordinates[0];
    double *p_density = &ws_density[0];
    double *p_vrtm = &ws_vrtm[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape function
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);
    
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      num_face_nodes = bulk_data.num_nodes(face);
      // sanity check on num nodes
      ThrowAssert( num_face_nodes == nodesPerFace );
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        // gather scalars
        p_density[ni] = *stk::mesh::field_data(*density_, node);
        // gather vectors
        const double * vrtm = stk::mesh::field_data(*velocityRTM_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_vrtm[offSet+j] = vrtm[j];
        }
      }
      
      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const * face_elem_rels = b.begin_elements(k);
      ThrowAssert( b.num_elements(k) == 1 );

      // get element; its face ordinal number and populate face_node_ordinals
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = b.begin_element_ordinals(k)[0];
      theElemTopo.side_node_ordinals(face_ordinal, face_node_ordinals.begin());

      //==========================================
      // gather nodal data off of element
      //==========================================
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);
      int num_nodes = bulk_data.num_nodes(element);
      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        // gather vectors
        double * coords = stk::mesh::field_data(*coordinates_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
        }
      }

      for ( int ip = 0; ip < num_face_nodes; ++ip ) {

        const int opposingNode = meSCS->opposingNodes(face_ordinal,ip);
        const int nearestNode = face_node_ordinals[ip];

        // node on the face
        stk::mesh::Entity nodeR = elem_node_rels[nearestNode];

        // pointers to nearest node data
        double *ncNormalFlux = stk::mesh::field_data(*ncNormalFlux_, nodeR);
        double *ncPenalty = stk::mesh::field_data(*ncPenalty_, nodeR);
        double *ncArea = stk::mesh::field_data(*ncArea_, nodeR);
       
        // offset for bip area vector and types of shape function
        const int faceOffSet = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        // interpolate to bip; sneak in aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          p_uBip[j] = 0.0;
          p_rho_uBip[j] = 0.0;
          const double axj = areaVec[faceOffSet+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // for bip values
        double rhoBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];    
          const double rhoIC = p_density[ic];
          rhoBip += r*rhoIC;
          const int offSetFN = ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_uBip[j] += r*p_vrtm[offSetFN+j];
            p_rho_uBip[j] += r*rhoIC*p_vrtm[offSetFN+j];
          }
        }        

        double mdot = 0.0;
        for ( int j = 0; j < nDim; ++j )
          mdot += (interpTogether*p_rho_uBip[j] + om_interpTogether*rhoBip*p_uBip[j])*areaVec[faceOffSet+j];

        // characteristic length
        double charLength = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          double dxj = p_coordinates[opposingNode*nDim+j] - p_coordinates[nearestNode*nDim+j];
          charLength += dxj*dxj;
        }
        charLength = std::sqrt(charLength);

        // assemble the nodal quantities
        *ncNormalFlux += mdot;
        *ncPenalty += projTimeScale/charLength*aMag;
        *ncArea += aMag;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleContinuityNonConformalPenaltyAlgorithm::~AssembleContinuityNonConformalPenaltyAlgorithm()
{
  // does nothing
}

} // namespace nalu
} // namespace Sierra
