/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef HaloInfo_h
#define HaloInfo_h

//==============================================================================
// Includes and forwards
//==============================================================================

#include <stk_mesh/base/Entity.hpp>

#include <vector>

namespace sierra {
namespace nalu {

//=============================================================================
// Class Definition
//=============================================================================
// HaloInfo
//=============================================================================
/**
 * * @par Description:
 * - class for halo stuff.
 *
 * @par Design Considerations:
 * -
 */
//=============================================================================
class HaloInfo {

 public:

  // constructor and destructor
  HaloInfo(
    stk::mesh::Entity node,
    const int nDim );

  ~HaloInfo();

  void resize_vectors(const int nDim);
  
  stk::mesh::Entity faceNode_;
  stk::mesh::Entity owningElement_;
  stk::mesh::Entity prevOwningElement_;
  
  std::vector<double> haloEdgeAreaVec_;
  std::vector<double> haloNodalCoords_;
  std::vector<double> haloMeshVelocity_;
  std::vector<double> checkhaloNodalCoords_;
  std::vector<double> nodalCoords_;
  std::vector<double> isoParCoords_;
  double haloEdgeDs_;
  double bestX_;
  int elemIsGhosted_;
  
};
  
} // end sierra namespace
} // end nalu namespace

#endif
