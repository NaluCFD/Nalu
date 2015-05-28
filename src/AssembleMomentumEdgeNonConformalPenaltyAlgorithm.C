/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleMomentumEdgeNonConformalPenaltyAlgorithm.h>

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
// AssembleMomentumEdgeNonConformalPenaltyAlgorithm - nodal lambda, flux
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumEdgeNonConformalPenaltyAlgorithm::AssembleMomentumEdgeNonConformalPenaltyAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  VectorFieldType *velocity,
  VectorFieldType *ncNormalFlux,
  ScalarFieldType *ncPenalty,
  ScalarFieldType *ncArea,
  ScalarFieldType *diffFluxCoeff)
  : Algorithm(realm, part),
    includeDivU_(realm_.get_divU()),
    velocity_(velocity),
    ncNormalFlux_(ncNormalFlux),
    ncPenalty_(ncPenalty),
    ncArea_(ncArea),
    diffFluxCoeff_(diffFluxCoeff),
    dudx_(NULL),
    coordinates_(NULL),
    exposedAreaVec_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeNonConformalPenaltyAlgorithm::execute()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // space for dui/dxj; the modified gradient with NOC
  std::vector<double> duidxj(nDim*nDim);
  std::vector<double> fx(nDim);

  // pointers to fixed values
  double *p_duidxj = &duidxj[0];
  double *p_fx = &fx[0];

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);

  // nodal fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_velocity;
  std::vector<double> ws_diffFluxCoeff;
 
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
    const int numScsIps = meFC->numIntPoints_;

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

      // get the relations; populate connected nodes
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);

      for ( int ip = 0; ip < numScsIps; ++ip ) {

        const int opposingNode = meSCS->opposingNodes(face_ordinal,ip);
        const int nearestNode = face_node_ordinals[ip];

        // left and right nodes; right is on the face; left is the opposing node
        stk::mesh::Entity nodeL = elem_node_rels[opposingNode];
        stk::mesh::Entity nodeR = elem_node_rels[nearestNode];

        // extract nodal fields
        const double * coordL = stk::mesh::field_data(*coordinates_, nodeL);
        const double * coordR = stk::mesh::field_data(*coordinates_, nodeR);

        const double * uNp1L = stk::mesh::field_data(velocityNp1, nodeL);
        const double * uNp1R = stk::mesh::field_data(velocityNp1, nodeR);

        const double diffFluxCoeffBip = *stk::mesh::field_data(*diffFluxCoeff_, nodeR);
        const double *dudxR = stk::mesh::field_data(*dudx_, nodeR);

        // pointers to nearest node data
        double *ncNormalFlux = stk::mesh::field_data(*ncNormalFlux_, nodeR);
        double *ncPenalty = stk::mesh::field_data(*ncPenalty_, nodeR);
        double *ncArea = stk::mesh::field_data(*ncArea_, nodeR);
       
        // offset for bip area vector and types of shape function
        const int faceOffSet = ip*nDim;

        // compute geometry and characteristic length
        double axdx = 0.0;
        double asq = 0.0;
        double charLength = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[faceOffSet+j];
          const double dxj = coordR[j]  - coordL[j];
          asq += axj*axj;
          axdx += axj*dxj;
          charLength += dxj*dxj;
        }
        const double inv_axdx = 1.0/axdx;
        const double aMag = std::sqrt(asq);
        charLength = std::sqrt(charLength);

        // form duidxj with over-relaxed procedure of Jasak:
        for ( int i = 0; i < nDim; ++i ) {

          // difference between R and L nodes for component i
          const double uidiff = uNp1R[i] - uNp1L[i];

          // offset into all forms of dudx
          const int offSetI = nDim*i;

          // start sum for NOC contribution
          double GlUidxl = 0.0;
          for ( int l = 0; l< nDim; ++l ) {
            const int offSetIL = offSetI+l;
            const double dxl = coordR[l] - coordL[l];
            const double GlUi = dudxR[offSetIL];
            GlUidxl += GlUi*dxl;
          }

          // form full tensor dui/dxj with NOC
          for ( int j = 0; j < nDim; ++j ) {
            const int offSetIJ = offSetI+j;
            const double axj = areaVec[faceOffSet+j];
            const double GjUi = dudxR[offSetIJ];
            p_duidxj[offSetIJ] = GjUi + (uidiff - GlUidxl)*axj*inv_axdx;
          }
        }

        // divU
        double divU = 0.0;
        for ( int j = 0; j < nDim; ++j)
          divU += p_duidxj[j*nDim+j];

        // form viscous stress
        for (int i = 0; i < nDim; ++i ) {
          const double axi = areaVec[faceOffSet+i];
          double fxi = 2.0/3.0*diffFluxCoeffBip*divU*axi*includeDivU_;
          const int offSetI = nDim*i;
          for ( int j = 0; j < nDim; ++j ) {
            const int offSetTrans = nDim*j+i;
            const double axj = areaVec[faceOffSet+j];
            fxi += -diffFluxCoeffBip*(p_duidxj[offSetI+j] + p_duidxj[offSetTrans])*axj;
          }
          p_fx[i] = fxi;
        }

        // assemble the nodal quantities; scalar
        *ncArea += aMag;
        *ncPenalty += diffFluxCoeffBip/charLength*aMag;

        // vector
        for ( int i = 0; i < nDim; ++i ) {
          ncNormalFlux[i] += p_fx[i];
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumEdgeNonConformalPenaltyAlgorithm::~AssembleMomentumEdgeNonConformalPenaltyAlgorithm()
{
  // does nothing
}

} // namespace nalu
} // namespace Sierra
