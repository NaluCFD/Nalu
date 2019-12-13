/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeWallFrictionVelocityProjectedAlgorithm_h
#define ComputeWallFrictionVelocityProjectedAlgorithm_h

#include<Algorithm.h>
#include<PointInfo.h>
#include<FieldTypeDef.h>

// stk
#include <stk_search/SearchMethod.hpp>

// c++
#include <vector>

// forwards
namespace stk {
  namespace mesh {
    class BulkData;
    class Ghosting;
    class MetaData;
    class Part;
  }
}

namespace sierra{
namespace nalu{

class Realm;

class ComputeWallFrictionVelocityProjectedAlgorithm : public Algorithm
{
public:

  ComputeWallFrictionVelocityProjectedAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const double projectedDistance,
    const bool useShifted,
    std::vector<std::vector<PointInfo *> > &pointInfoVec,
    stk::mesh::Ghosting *wallFunctionGhosting);
  virtual ~ComputeWallFrictionVelocityProjectedAlgorithm();

  void execute();

  void set_data( 
    double theDouble);

  void compute_utau(
      const double &up, const double &yp,
      const double &density, const double &viscosity,
      double &utau);
  
  // debug
  void provide_output(const PointInfo *pInfo, const bool problemPoint);
  
  // ghosting procedure set of calls
  void initialize();
  void initialize_ghosting();
  void construct_bounding_points();
  void reset_point_info();
  void construct_bounding_boxes();
  void coarse_search();
  void manage_ghosting();
  void complete_search();
  
  const bool useShifted_;
  /* vector of PointInfo pointInfoVec_[k] provides numScsBip PointInfo objects on face k */
  std::vector<std::vector<PointInfo *> > &pointInfoVec_;  
  stk::mesh::Ghosting *wallFunctionGhosting_;

  stk::mesh::BulkData *bulkData_;
  stk::mesh::MetaData *metaData_;
  const int nDim_;

  const double yplusCrit_;
  const double elog_;
  const double kappa_;
  const int maxIteration_;
  const double tolerance_;

  bool firstInitialization_;
  const bool provideOutput_;
 
  /* search data structures that need to be initialized*/
  const stk::search::SearchMethod searchMethod_;
  const double expandBoxPercentage_;
  size_t needToGhostCount_;

  VectorFieldType *velocity_;
  VectorFieldType *bcVelocity_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallNormalDistanceBip_;
  ScalarFieldType *assembledWallNormalDistance_;
  ScalarFieldType *assembledWallArea_;

  // per part value for projected distance
  std::vector<double> projectedDistanceVec_;

  /* data types for stk_search */
  std::vector<boundingPoint> boundingPointVec_;
  std::vector<boundingBox> boundingBoxVec_;
  std::vector<std::pair<uint64IdentProc, uint64IdentProc> > searchKeyPair_;

  // vector of elements to ghost
  stk::mesh::EntityProcVec elemsToGhost_;
  
  // data structure to parallel communicate nodal data to ghosted elements
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
