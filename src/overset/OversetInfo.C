/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <overset/OversetInfo.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// OversetInfo - contains orphan point -> donor elements
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
OversetInfo::OversetInfo(
  stk::mesh::Entity node,
  const int nDim)
  : orphanNode_(node),
    owningElement_(),
    bestX_(1.0e16),
    elemIsGhosted_(0),
    meSCS_(NULL)
{
  // resize stuff
  isoParCoords_.resize(nDim);
  nodalCoords_.resize(nDim);
}
//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
OversetInfo::~OversetInfo()
{
  // nothing to delete
}

} // namespace NaluUnit
} // namespace sierra
