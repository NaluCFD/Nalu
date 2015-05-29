/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleScalarEdgeNonConformalPenaltyAlgorithm.h>

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
// AssembleScalarEdgeNonConformalPenaltyAlgorithm - nodal lambda, diff flux
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEdgeNonConformalPenaltyAlgorithm::AssembleScalarEdgeNonConformalPenaltyAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *scalarQ,
  VectorFieldType *GjQ,
  ScalarFieldType *ncNormalFlux,
  ScalarFieldType *ncPenalty,
  ScalarFieldType *ncArea,
  ScalarFieldType *diffFluxCoeff)
  : Algorithm(realm, part),
    scalarQ_(scalarQ),
    GjQ_(GjQ),
    ncNormalFlux_(ncNormalFlux),
    ncPenalty_(ncPenalty),
    ncArea_(ncArea),
    diffFluxCoeff_(diffFluxCoeff),
    coordinates_(NULL),
    exposedAreaVec_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEdgeNonConformalPenaltyAlgorithm::execute()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

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
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);
     
      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const * face_elem_rels = b.begin_elements(k);
      ThrowAssert( b.num_elements(k) == 1 );

      // get element; its face ordinal number and populate face_node_ordinals
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = b.begin_element_ordinals(k)[0];
      theElemTopo.side_node_ordinals(face_ordinal, face_node_ordinals.begin());
      
      // get the relations
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);

      for ( int ip = 0; ip < num_face_nodes; ++ip ) {

        const int opposingNode = meSCS->opposingNodes(face_ordinal,ip);
        const int nearestNode = face_node_ordinals[ip];

        // left and right nodes; right is on the face; left is the opposing node
        stk::mesh::Entity nodeL = elem_node_rels[opposingNode];
        stk::mesh::Entity nodeR = elem_node_rels[nearestNode];

        // extract nodal fields
        const double * coordL = stk::mesh::field_data(*coordinates_, nodeL );
        const double * coordR = stk::mesh::field_data(*coordinates_, nodeR );

        const double scalarQL = *stk::mesh::field_data(*scalarQ_, nodeL );
        const double scalarQR = *stk::mesh::field_data(*scalarQ_, nodeR );

        // pointers to nearest node data; const and data to be assembled
        const double * GjQR    =  stk::mesh::field_data(*GjQ_, nodeR );
        const double diffFluxCoeffR = *stk::mesh::field_data(*diffFluxCoeff_, nodeR );
        double *ncNormalFlux = stk::mesh::field_data(*ncNormalFlux_, nodeR);
        double *ncPenalty = stk::mesh::field_data(*ncPenalty_, nodeR);
        double *ncArea = stk::mesh::field_data(*ncArea_, nodeR);
       
        // offset for bip area vector and types of shape function
        const int faceOffSet = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        // compute geometry
        double axdx = 0.0;
        double asq = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[faceOffSet+j];
          const double dxj = coordR[j] - coordL[j];
          asq += axj*axj;
          axdx += axj*dxj;
        }
        const double inv_axdx = 1.0/axdx;
        const double aMag = std::sqrt(asq);

        // NOC and characteristic length
        double nonOrth = 0.0;
        double charLength = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[faceOffSet+j];
          const double dxj = coordR[j] - coordL[j];
          const double kxj = axj - asq*inv_axdx*dxj;
          nonOrth += -diffFluxCoeffR*kxj*GjQR[j];
          charLength += dxj*dxj;
        }
        charLength = std::sqrt(charLength);
                
        // assemble the nodal quantities
        *ncNormalFlux += -diffFluxCoeffR*(scalarQR-scalarQL)*asq*inv_axdx + nonOrth;
        *ncPenalty += diffFluxCoeffR/charLength*aMag;
        *ncArea += aMag;
      }
    }
  }  
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEdgeNonConformalPenaltyAlgorithm::~AssembleScalarEdgeNonConformalPenaltyAlgorithm()
{
  // does nothing
}

} // namespace nalu
} // namespace Sierra
