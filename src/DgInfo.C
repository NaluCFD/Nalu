/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <DgInfo.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Acon_DgInfo - contains non-conformal DG-based information
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
DgInfo::DgInfo(
  int parallelRank,
  uint64_t globalFaceId,
  uint64_t localGaussPointId,
  int currentGaussPointId,
  stk::mesh::Entity currentFace,
  stk::mesh::Entity currentElement,
  const int currentFaceOrdinal,
  MasterElement *meFCCurrent,
  MasterElement *meSCSCurrent,
  stk::topology currentElementTopo,
  const int nDim)
  : parallelRank_(parallelRank),
    globalFaceId_(globalFaceId),
    localGaussPointId_(localGaussPointId),
    currentGaussPointId_(currentGaussPointId),
    currentFace_(currentFace),
    currentElement_(currentElement),
    currentFaceOrdinal_(currentFaceOrdinal),
    meFCCurrent_(meFCCurrent),
    meSCSCurrent_(meSCSCurrent),
    currentElementTopo_(currentElementTopo),
    nDim_(nDim),
    bestX_(1.0e16),
    opposingFaceIsGhosted_(0)
{
  // resize internal vectors
  currentGaussPointCoords_.resize(nDim);
  // isoPar coords will map to full volume element
  currentIsoParCoords_.resize(nDim);
  opposingIsoParCoords_.resize(nDim);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
DgInfo::~DgInfo()
{
  // nothing to delete
}

} // namespace Acon
} // namespace sierra
