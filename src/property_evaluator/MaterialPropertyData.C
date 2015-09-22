/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/MaterialPropertyData.h>

#include <Enums.h>

namespace sierra{
namespace nalu{

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MaterialPropertyData::MaterialPropertyData()
  : type_(MaterialPropertyType_END),
    constValue_(0.0),
    primary_(0.0),
    secondary_(0.0),
    auxVarName_("na"),
    tablePropName_("na"),
    tableAuxVarName_("na"),
    genericPropertyEvaluatorName_("na")
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MaterialPropertyData::~MaterialPropertyData()
{
  // nothing
}

} // namespace nalu
} // namespace Sierra
