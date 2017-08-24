/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MasterElementFactory_h
#define MasterElementFactory_h

#include <string>

namespace stk { struct topology; }

namespace sierra{
namespace nalu{

  class MasterElement;

  MasterElement*
  get_surface_master_element(
    const stk::topology& theTopo,
    int dimension = 0,
    std::string quadType = "GaussLegendre");

  MasterElement*
  get_volume_master_element(
    const stk::topology& theTopo,
    int dimension = 0,
    std::string quadType = "GaussLegendre");
    
} // namespace nalu
} // namespace Sierra

#endif
