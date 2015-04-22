/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNodalGradUElemContactAlgorithm.h>
#include <SolverAlgorithm.h>

#include <ContactInfo.h>
#include <ContactManager.h>
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
// AssembleNodalGradUElemContactAlgorithm - contact scalar dqdx
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradUElemContactAlgorithm::AssembleNodalGradUElemContactAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *vectorQ,
    GenericFieldType *dqdx,
    VectorFieldType *haloQ)
  : Algorithm(realm, part),
    vectorQ_(vectorQ),
    dqdx_(dqdx),
    haloQ_(haloQ)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // populate fieldVec
  ghostFieldVec_.push_back(vectorQ_);

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradUElemContactAlgorithm::execute()
{
  
  // determine halo state; save this value to the face node
  populate_halo_state();

  // now add in elemental projected nodal gradient contributions
  add_elem_gradq();
  
}

//--------------------------------------------------------------------------
//-------- populate_halo_state ---------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradUElemContactAlgorithm::populate_halo_state()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // parallel communicate ghosted entities
  if ( NULL != realm_.contactManager_->contactGhosting_ )
    stk::mesh::communicate_field_data(*(realm_.contactManager_->contactGhosting_), ghostFieldVec_);
  
  // sizing
  const int nDim = meta_data.spatial_dimension();
  std::vector<double> vectorQR(nDim);

  // iterate contactInfoVec_
  std::vector<ContactInfo *>::iterator ii;
  for( ii=realm_.contactManager_->contactInfoVec_.begin();
       ii!=realm_.contactManager_->contactInfoVec_.end(); ++ii ) {

    // get master element type for this contactInfo
    MasterElement *meSCS  = (*ii)->meSCS_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    std::vector <double > elemNodalVectorQ(nodesPerElement*nDim);
    std::vector <double > shpfc(nodesPerElement);

    // iterate halo face nodes
    std::map<uint64_t, HaloInfo *>::iterator iterHalo;
    for (iterHalo  = (*ii)->haloInfoMap_.begin();
         iterHalo != (*ii)->haloInfoMap_.end();
         ++iterHalo) {

      // halo info object of interest
      HaloInfo * infoObject = (*iterHalo).second;

      // extract element mesh object and global id for face node
      stk::mesh::Entity elem  = infoObject->owningElement_;

      stk::mesh::Entity const *elem_node_rels = bulk_data.begin_nodes(elem);
      const unsigned num_nodes = bulk_data.num_nodes(elem);
      
      // now load the elemental values for future interpolation
      for ( unsigned ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        const double *vectorQ = stk::mesh::field_data(*vectorQ_, node);
        for ( int j = 0; j < nDim; ++j )
          elemNodalVectorQ[ni*nDim+j] = vectorQ[j];
      }

      // extract nodal fields; right state is Halo and requires inperpolation
      const int sizeOfScalarField = nDim;
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalVectorQ[0],
        &vectorQR[0]);
      
      double *haloQ  = stk::mesh::field_data(*haloQ_, infoObject->faceNode_);
      for (int j = 0; j < nDim; ++j)
        haloQ[j] = vectorQR[j];
      
    }
  }

  // parallel communicate locally owned nodes to shared
  std::vector< const stk::mesh::FieldBase *> fieldVec(1, haloQ_);
  stk::mesh::copy_owned_to_shared( bulk_data, fieldVec);
  
}

//--------------------------------------------------------------------------
//-------- add_elem_gradq --------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradUElemContactAlgorithm::add_elem_gradq()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  VectorFieldType *haloDxj = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_dxj");
 
  const int nDim = meta_data.spatial_dimension();

  // loop over locally owned faces and construct missing elemental contributions
  stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);
  
  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned );  
  
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
	ib != face_buckets.end() ; ++ib ) {
    
    stk::mesh::Bucket & b = **ib;
    
    // extract master element; hard coded for quad or hex; 
    // quad is always true for 2D while for 3D, either hex or wedge apply
    const stk::topology & theElemTopo = (nDim == 2) ? stk::topology::QUAD_4_2D : stk::topology::HEX_8;
    const int num_face_nodes = (nDim == 2) ? 2 : 4;
    std::vector<int> face_node_ordinals(num_face_nodes);
    
    // extract master element for extruded element type
    MasterElement *meSCS = realm_.get_surface_master_element(theElemTopo);
    MasterElement *meSCV = realm_.get_volume_master_element(theElemTopo);
    
    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;
    
    // mapping between exposed face and extruded element's overlapping face
    const int *faceNodeOnExtrudedElem = meSCS->faceNodeOnExtrudedElem();

    // mapping between exposed face and extruded element's opposing face
    const int *opposingNodeOnExtrudedElem = meSCS->opposingNodeOnExtrudedElem();
  
    // mapping between exposed face scs ips and halo edge
    const int *faceScsIpOnExtrudedElem = meSCS->faceScsIpOnExtrudedElem();
  
    // mapping between exposed face scs ips and exposed face edge
    const int *faceScsIpOnFaceEdges = meSCS->faceScsIpOnFaceEdges();
    
    // alignment of face:edge ordering and scsip area vector
    const double *edgeAlignedArea = meSCS->edgeAlignedArea();

    // define scratch field
    std::vector<double > ws_coordinates(nodesPerElement*nDim);
    std::vector<double > ws_scs_areav(numScsIp*nDim);
    std::vector<double > ws_vectorQ(nodesPerElement*nDim);
    std::vector<double> ws_shape_function(numScsIp*nodesPerElement);

    // pointers
    double *p_shape_function = &ws_shape_function[0];
    meSCS->shape_fcn(&p_shape_function[0]);

    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get face
      stk::mesh::Entity face = b[k];
      
      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const* face_elem_rels = bulk_data.begin_elements(face);
      stk::mesh::ConnectivityOrdinal const* face_elem_ords = bulk_data.begin_element_ordinals(face);
      const int num_elements = bulk_data.num_elements(face);
      ThrowRequire( num_elements == 1 );
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = face_elem_ords[0];
      theElemTopo.side_node_ordinals(face_ordinal, face_node_ordinals.begin());
      
      // concentrate on loading up the nodal coordinates/vectorQ for the extruded element
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        const double * coords = stk::mesh::field_data(*coordinates, node);
        const double * hDxj = stk::mesh::field_data( *haloDxj, node );
        const double * haloQ = stk::mesh::field_data( *haloQ_, node );
        const double * vectorQ = stk::mesh::field_data( *vectorQ_, node );
        const int faceNode = faceNodeOnExtrudedElem[face_ordinal*num_nodes + ni];
        const int opposingNode = opposingNodeOnExtrudedElem[face_ordinal*num_nodes + ni];
        const int offSetFN = faceNode*nDim;
        const int offSetON = opposingNode*nDim;

        // populate vectors; face node and halo node
        for ( int j=0; j < nDim; ++j ) {
          ws_coordinates[offSetFN+j] = coords[j];
          ws_vectorQ[offSetFN+j] = vectorQ[j];
          ws_coordinates[offSetON+j] = coords[j] + hDxj[j];
          ws_vectorQ[offSetON+j] = haloQ[j];
        }
      }      
      
      // compute scs integration point areavec
      double scs_error = 0.0;
      meSCS->determinant(1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);
      
      // assemble halo ip contribution for face node
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        const double &dualNodalVolume = *stk::mesh::field_data( *dualNodalVolume_, node );
        // area vector for halo edge;
        // face ordinal 0 for extruded element has all scs area vectors pointing from face to opposing face
        const int scsIp = faceScsIpOnExtrudedElem[face_ordinal*num_nodes + ni];

        // interpolate element nodal values to this scsIp of interest
        double qIp[3] = {0.0,0.0,0.0};
        const int offSetIp = scsIp*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function[offSetIp + ic];
          for (int j = 0; j < nDim; ++j ) {
            qIp[j] += r*ws_vectorQ[ic*nDim+j];
          }
        }
	
        // add in nodal gradient contribution
        double *dqdx = stk::mesh::field_data( *dqdx_, node );
        for ( int i = 0; i < nDim; ++i ) {
          const int row_gradQ = i*nDim;
          const double qip = qIp[i];
          for ( int j = 0; j < nDim; ++j ) {
            dqdx[row_gradQ+j] += qip*ws_scs_areav[scsIp*nDim+j]/dualNodalVolume;
          }
        }
      }
      
      // deal with edges on the exposed face and each
      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(element);
      
      // face edge relations; if this is 2D then the face is a edge and size is unity
      stk::mesh::Entity const* face_edge_rels = bulk_data.begin_edges(face);
      const int num_face_edges = bulk_data.num_edges(face);
      
      int num_edges = (nDim == 3) ? num_face_edges : 1;
      
      for ( int ie = 0; ie < num_edges; ++ie ) {
	
        // get edge
        stk::mesh::Entity edge = (nDim == 3) ? face_edge_rels[ie] : face;
	
        // get the relations from edge
        stk::mesh::Entity const* edge_node_rels = bulk_data.begin_nodes(edge);
        const int edge_num_nodes = bulk_data.num_nodes(edge);
        // sanity check on num nodes
        if ( edge_num_nodes != 2 ){
          throw std::runtime_error("num nodes is not 2");
        }
	
        // extract ip for this edge
        const int scsIp = faceScsIpOnFaceEdges[face_ordinal*num_edges + ie];
	
        // correct area for edge and scs area vector from extruded element alignment
        const double alignmentFac = edgeAlignedArea[face_ordinal*num_edges + ie];
	
        // interpolate element nodal values to this scsIp of interest
        double qIp[3] = {0.0,0.0,0.0};
        const int offSetIp = scsIp*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function[offSetIp + ic];
          for (int j = 0; j < nDim; ++j ) {
            qIp[j] += r*ws_vectorQ[ic*nDim+j];
          }
        }
	
        // left and right nodes on the edge
        stk::mesh::Entity nodeL = edge_node_rels[0];
        stk::mesh::Entity nodeR = edge_node_rels[1];
	
        // does edge point correctly
        const int leftNode = face_node_ordinals[ie];
        const size_t iglob_Lelem = bulk_data.identifier(elem_node_rels[leftNode]);
        const size_t iglob_Ledge = bulk_data.identifier(edge_node_rels[0]);
	
        // determine the sign value for area vector; if Left node is the same,
        // then the element and edge relations are aligned
        const double sign = ( iglob_Lelem == iglob_Ledge ) ? 1.0 : -1.0;
	
        // add in nodal gradient contribution
        double *dqdxL = stk::mesh::field_data( *dqdx_, nodeL );
        double *dqdxR = stk::mesh::field_data( *dqdx_, nodeR );
        const double &dualNodalVolumeL = *stk::mesh::field_data( *dualNodalVolume_, nodeL );
        const double &dualNodalVolumeR = *stk::mesh::field_data( *dualNodalVolume_, nodeR );
        for (int i = 0; i < nDim; ++i ) {
          const int row_gradQ = i*nDim;
          const double qip = qIp[i];
          for ( int j = 0; j < nDim; ++j ) {
            dqdxL[row_gradQ+j] += qip*ws_scs_areav[scsIp*nDim+j]/dualNodalVolumeL*sign*alignmentFac;
            dqdxR[row_gradQ+j] -= qip*ws_scs_areav[scsIp*nDim+j]/dualNodalVolumeR*sign*alignmentFac;
          }
	  
        }
      }
    }
  }

  // parallel assembly handled elsewhere
  
}

} // namespace nalu
} // namespace Sierra
