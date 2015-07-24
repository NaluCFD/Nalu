/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef OversetInfo_h
#define OversetInfo_h

//==============================================================================
// Includes and forwards
//==============================================================================

#include <stk_mesh/base/Entity.hpp>
#include <cmath> 
#include <vector>

namespace sierra {
namespace nalu {

class MasterElement;

//=============================================================================
// Class Definition
//=============================================================================
// OversetInfo
//=============================================================================
class OversetInfo {

 public:

  // constructor and destructor
  OversetInfo(
    stk::mesh::Entity node,
    const int nDim );

  ~OversetInfo();

  stk::mesh::Entity orphanNode_;
  stk::mesh::Entity owningElement_;

  double bestX_;
  int elemIsGhosted_;

  // master element for background mesh
  MasterElement *meSCS_;

  std::vector<double> isoParCoords_;
  std::vector<double> nodalCoords_;

};
  
} // end sierra namespace
} // end naluUnit namespace

#endif
