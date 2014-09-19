/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LinInterp_h
#define LinInterp_h


#include <string>
#include <vector>
#include <utility>

#include <stk_mesh/base/Entity.hpp>

#include <Realm.h>
#include <master_element/MasterElement.h>

namespace sierra{
namespace nalu{

template <class FROM, class TO> class LinInterp {
  
public :

  typedef FROM                            MeshA;
  typedef TO                              MeshB;
  typedef typename MeshA::EntityKey       EntityKeyA;
  typedef typename MeshB::EntityKey       EntityKeyB;
  typedef typename MeshA::EntityProc      EntityProcA;
  typedef typename MeshB::EntityProc      EntityProcB;
  
  typedef std::pair<EntityProcB, EntityProcA>          EntityProcRelation;
  typedef std::vector<EntityProcRelation>              EntityProcRelationVec;

  typedef std::multimap<EntityKeyB, EntityKeyA>        EntityKeyMap;
  
  enum { Dimension = MeshA::Dimension };

  static void filter_to_nearest(EntityKeyMap    &RangeToDomain,
      const MeshA     &FromElem,
      MeshB     &ToPoints) ;
    
  static void apply (MeshB         &ToPoints,
      const MeshA         &FromElem,
      const EntityKeyMap &RangeToDomain) ;
};

template <class FROM, class TO>  void LinInterp<FROM,TO>::filter_to_nearest (
  EntityKeyMap    &RangeToDomain,
  const MeshA     &FromElem,
  MeshB           &ToPoints) {

  const stk::mesh::BulkData &fromBulkData = FromElem.fromBulkData_;
  stk::mesh::BulkData &toBulkData = ToPoints.toBulkData_;
  Realm &fromRealm = FromElem.fromRealm_;

  const VectorFieldType *fromcoordinates = FromElem.fromcoordinates_;
  const VectorFieldType *tocoordinates   = ToPoints.tocoordinates_;

  // FixMe: Check nDim against Dimension to make sure they are equal
  const unsigned nDim = FromElem.fromMetaData_.spatial_dimension();

  typedef typename EntityKeyMap::iterator iterator;
  typedef typename EntityKeyMap::const_iterator const_iterator;

  for (const_iterator current_key=RangeToDomain.begin(); current_key!=RangeToDomain.end(); ) { 

    double bestX_ = std::numeric_limits<double>::max();

    const stk::mesh::EntityKey thePt  = current_key->first;
    stk::mesh::Entity theNode = toBulkData.get_entity(thePt);
    // load nodal coordinates from node
    const double * tocoords = stk::mesh::field_data(*tocoordinates, theNode );

    std::pair<iterator, iterator> keys=RangeToDomain.equal_range(current_key->first);
    iterator nearest = keys.second;

    for (iterator ii=keys.first; ii != keys.second; ++ii) {

      const stk::mesh::EntityKey theBox = ii->second; 
      stk::mesh::Entity theElem = fromBulkData.get_entity(theBox);    
    
      // extract master element from the bucket in which the element resides
      const stk::mesh::Bucket &theBucket = fromBulkData.bucket(theElem);
      const stk::topology &theElemTopo = theBucket.topology();
      MasterElement *meSCS = fromRealm.get_surface_master_element(theElemTopo);

      // load nodal coordinates from element
      stk::mesh::Entity const* elem_node_rels = fromBulkData.begin_nodes(theElem);
      const int num_nodes = fromBulkData.num_nodes(theElem);
    
      const int nodesPerElement = meSCS->nodesPerElement_;
      std::vector<double> theElementCoords(nDim*nodesPerElement);

      for ( int ni = 0; ni < num_nodes; ++ni ) { 
        stk::mesh::Entity node = elem_node_rels[ni];
     
        // load up vectors
        const double * fromcoords = stk::mesh::field_data(*fromcoordinates, node );
        for ( unsigned j = 0; j < nDim; ++j ) { 
          const int offSet = j*nodesPerElement + ni; 
          theElementCoords[offSet] = fromcoords[j];
        }   
      }   
    
      std::vector<double> isoParCoords(nDim);
      const double nearestDistance = meSCS->isInElement(&theElementCoords[0],
                                                        &(tocoords[0]),
                                                        &(isoParCoords[0]));
      if ( nearestDistance < bestX_ ) { 
        bestX_         = nearestDistance;    
        ToPoints.TransferInfo_[thePt] = isoParCoords;
        nearest = ii;
      }   
    }
    current_key = keys.second;
    if (nearest != keys.first ) RangeToDomain.erase(keys.first, nearest);
    if (nearest != keys.second) RangeToDomain.erase(++nearest, keys.second);
  }
}


template <class FROM, class TO>  void LinInterp<FROM,TO>::apply 
       (MeshB              &ToPoints,
        const MeshA        &FromElem,
        const EntityKeyMap &RangeToDomain) {
  
  const stk::mesh::BulkData &fromBulkData = FromElem.fromBulkData_;
  stk::mesh::BulkData         &toBulkData = ToPoints.toBulkData_;
  Realm &fromRealm = FromElem.fromRealm_;

  const int sizeOfScalarField = 1;

  typename EntityKeyMap::const_iterator ii;
  for(ii=RangeToDomain.begin(); ii!=RangeToDomain.end(); ++ii ) { 
    
    const stk::mesh::EntityKey thePt  = ii->first;
    const stk::mesh::EntityKey theBox = ii->second; 
    
    if (1 != ToPoints.TransferInfo_.count(thePt)) {
      if (0 == ToPoints.TransferInfo_.count(thePt)) 
        throw std::runtime_error("Key not found in database");
      else  
        throw std::runtime_error("Too many Keys found in database");
    }
    const std::vector<double> &isoParCoords_ = ToPoints.TransferInfo_[thePt];
    stk::mesh::Entity theNode =   toBulkData.get_entity(thePt);
    stk::mesh::Entity theElem = fromBulkData.get_entity(theBox);    

    const stk::mesh::Bucket &theBucket = fromBulkData.bucket(theElem);
    const stk::topology &theElemTopo = theBucket.topology();
    MasterElement *meSCS = fromRealm.get_surface_master_element(theElemTopo);

    stk::mesh::Entity const* elem_node_rels = fromBulkData.begin_nodes(theElem);
    const int num_nodes = fromBulkData.num_nodes(theElem);
    const int nodesPerElement = meSCS->nodesPerElement_;

    std::vector <double > Coeff(nodesPerElement);

    for (unsigned n=0; n!=FromElem.fromFieldVec_.size(); ++n) {
      // now load the elemental values for future interpolation; fill in connected nodes
      for ( int ni = 0; ni < num_nodes; ++ni ) { 
        stk::mesh::Entity node = elem_node_rels[ni];
        Coeff[ni] = *stk::mesh::field_data(*FromElem.fromFieldVec_[n], node);
      }   

      double * toField = stk::mesh::field_data(*ToPoints.toFieldVec_[n], theNode);
      if (!toField) throw std::runtime_error("Receiving field undefined on mesh object.");
      meSCS->interpolatePoint(sizeOfScalarField,
                              &isoParCoords_[0],
                              &Coeff[0],
                              toField);
    }   
  }
}

} // namespace nalu
} // namespace Sierra

#endif
