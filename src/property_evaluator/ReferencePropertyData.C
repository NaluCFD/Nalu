/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/ReferencePropertyData.h>

#include <Enums.h>

namespace sierra{
namespace nalu{

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ReferencePropertyData::ReferencePropertyData()
  : speciesName_("na"),
    mw_(0.0),
    massFraction_(0.0),
    stoichiometry_(0.0),
    primaryMassFraction_(0.0),
    secondaryMassFraction_(0.0)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ReferencePropertyData::~ReferencePropertyData()
{
  // nothing
}

} // namespace nalu
} // namespace Sierra
