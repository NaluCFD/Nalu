/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <HaloInfo.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Acon_HaloInfo - contains virtual edge data
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HaloInfo::HaloInfo(
  stk::mesh::Entity node,
  const int nDim)
  : faceNode_(node),
    owningElement_(),
    prevOwningElement_(),
    bestX_(1.0e16),
    elemIsGhosted_(0)
{
  // resize stuff
  haloEdgeAreaVec_.resize(nDim);
  haloNodalCoords_.resize(nDim);
  haloMeshVelocity_.resize(nDim);
  checkhaloNodalCoords_.resize(nDim);
  nodalCoords_.resize(nDim);
  isoParCoords_.resize(nDim);
  
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
HaloInfo::~HaloInfo()
{
  // nothing to delete
}

} // namespace Acon
} // namespace sierra
