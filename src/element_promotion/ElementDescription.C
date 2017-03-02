/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level nalu      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <element_promotion/ElementDescription.h>

#include <element_promotion/QuadNElementDescription.h>
#include <element_promotion/HexNElementDescription.h>
#include <element_promotion/QuadratureRule.h>
#include <nalu_make_unique.h>

namespace sierra {
namespace nalu {

//--------------------------------------------------------------------------
std::unique_ptr<ElementDescription>
ElementDescription::create(int dimension, int order)
{
  auto nodeLocations1D = gauss_lobatto_legendre_rule(order + 1).first;
  if (dimension == 2) {
    return make_unique<QuadNElementDescription>(nodeLocations1D);
  }
  return make_unique<HexNElementDescription>(nodeLocations1D);
}
//--------------------------------------------------------------------------
ElementDescription::ElementDescription()  = default;
//--------------------------------------------------------------------------
ElementDescription::~ElementDescription() = default;
//--------------------------------------------------------------------------

} // namespace nalu
}  // namespace sierra
