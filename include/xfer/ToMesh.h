/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ToMesh_h
#define ToMesh_h

#include <string>
#include <vector>
#include <utility>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class ToMesh {
public :
  typedef stk::mesh:: Entity                             Entity;
  typedef std::vector<Entity>                            EntityVec;
  typedef stk::mesh:: EntityKey                          EntityKey;
  typedef std::set   <EntityKey>                         EntityKeySet;
  typedef stk::search::IdentProc<EntityKey, unsigned>    EntityProc;
  typedef std::vector<EntityProc>                        EntityProcVec;

  typedef stk::search::Point<double> Point;
  typedef stk::search::Sphere<double> Sphere;
  typedef std::pair<Sphere,EntityProc> BoundingBox;

  enum {Dimension = 3};

  typedef std::vector<std::pair<std::string, std::string> > PairNames;
  std::vector< const stk::mesh::FieldBase *>
  get_fields(const stk::mesh::MetaData  &toMetaData,
             const PairNames            &VarPairName) {
    std::vector< const stk::mesh::FieldBase *> toFieldVec;
    // provide field names
    for( PairNames::const_iterator i=VarPairName.begin(); i!=VarPairName.end(); ++i) {
      const std::string &name = i->second;
      const stk::mesh::FieldBase *tofield = stk::mesh::get_field_by_name(name,toMetaData);
      toFieldVec.push_back(tofield);
    }
    return toFieldVec;
  }

  ToMesh(stk::mesh::MetaData &toMetaData,
         stk::mesh::BulkData &toBulkData,
         const std::string &coordinates_name,
         const PairNames &VarPairName,
         const stk::mesh::Part  *toMeshPart,
         const stk::ParallelMachine comm,
         const double radius=.0001) :
    toMetaData_(toMetaData),
    toBulkData_(toBulkData),
    tocoordinates_(toMetaData.get_field<VectorFieldType>(stk::topology::NODE_RANK,coordinates_name)),
    toFieldVec_   (get_fields(toMetaData, VarPairName)),
    toMeshPart_(toMeshPart),
    comm_(comm),
    radius_(radius)   {}

  ~ToMesh(){};

  struct BoundingBoxCompare{
    bool operator()(const BoundingBox &a, const BoundingBox & b) const
    {
      return a.second.id() < b.second.id();
    }
  };


  // Needed for STK Transfer
  stk::ParallelMachine comm() const {return comm_;}

  void bounding_boxes (std::vector<BoundingBox> &v) const
  {

    const unsigned spatial_dimension = toMetaData_.spatial_dimension();

    Point center;

    stk::mesh::Selector s_locally_owned_union = toMetaData_.locally_owned_part()
      &stk::mesh::Selector(*toMeshPart_);

    stk::mesh::BucketVector const& node_buckets = toBulkData_.get_buckets( stk::topology::NODE_RANK, s_locally_owned_union );
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;

      const stk::mesh::Bucket::size_type length   = b.size();

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

        // get node
        stk::mesh::Entity node = b[k];

        double * coord = stk::mesh::field_data(*tocoordinates_, node);
        for (unsigned i=0; i<spatial_dimension; ++i) {
          center[i] = coord[i];
        }

        // setup ident
        ToMesh::EntityProc theIdent(toBulkData_.entity_key(node), toBulkData_.parallel_rank());

        ToMesh::BoundingBox theBox(Sphere(center,radius_), theIdent );
        v.push_back(theBox);
      }
    }
    std::sort(v.begin(), v.end(), BoundingBoxCompare());
  }

  void update_values()
  {
    std::vector<const stk::mesh::FieldBase *> fields(toFieldVec_.begin(), toFieldVec_.end());
    stk::mesh::copy_owned_to_shared  (  toBulkData_, fields);
  }

  stk::mesh::MetaData &toMetaData_;
  stk::mesh::BulkData &toBulkData_;
  const VectorFieldType     *tocoordinates_;
  const std::vector< const stk::mesh::FieldBase *> toFieldVec_;
  const stk::mesh::Part *toMeshPart_;
  const stk::ParallelMachine comm_;
  const double radius_;

  typedef std::map<stk::mesh::EntityKey, std::vector<double> > TransferInfo;
  TransferInfo TransferInfo_;

};

} // namespace nalu
} // namespace Sierra

#endif
